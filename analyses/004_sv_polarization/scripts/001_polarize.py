#!/usr/bin/env python3
"""Polarize SyRI SVs (indels >= 50 bp, INV/TRANS/INVTR) with outgroup genomes.

usage: 001_polarize.py <002 pairs.tsv> <002 syri dir> <cross syri dir> <out dir>

For every focal pair (F = ref, G = qry; rows of 002's pairs.tsv) each event is
expressed in F coordinates and checked against the F-vs-O alignment of every
outgroup O (O sets below).  Per outgroup the call is one of
  ins_F / del_F / ins_G / del_G   (indels: where the change happened)
  F / G                           (INV, TRANS, INVTR: lineage carrying the change)
  amb                             (not decidable from this outgroup)
and the consensus is the unique plurality of non-amb calls.

Outgroups: diploid x Cbp and diploid x diploid -> diploids of the *other* set
(Cr/Cg for Co pairs, Co for Cr pairs; the 9 cross pairs in 004 supply these);
Cbp x Cbp -> diploids of the same set (all exist in 002).
DUP / INVDP are counted but not polarized.

Writes <out dir>/01_polarized/<pair>.tsv and <out dir>/polarization_summary.tsv.
"""
import csv
import os
import sys
from bisect import bisect_right
from collections import Counter, defaultdict

pairs_tsv, dir002, dir_cross, out_dir = sys.argv[1:5]
MIN_SIZE = 50        # smallest event polarized
SMALL_INDEL = 20     # indels >= this are subtracted from "aligned" coverage
FRAC = 0.7           # coverage fraction to call has / lacks
ANCHOR_TOL = 50      # bp, O insertion vs G insertion position
SIZE_RATIO = 2.0     # O insertion must be within this factor of G's
FLANK = 100          # bp each side of an insertion point that must be aligned
SV_FRAC = 0.5        # fraction of F interval with O block type for INV/TRANS

DIP = {"Co": ["Co39", "CoMosher", "CoRussian"], "Cr": ["Cr145", "CrMosher", "CgMosher"]}
SET_OF = {g: s for s, gs in DIP.items() for g in gs}
CBP_SET = lambda g: g[-2:]   # CbpChineseCo -> Co
AL_TYPES = {"SYNAL", "INVAL", "TRANSAL", "INVTRAL", "DUPAL", "INVDPAL"}   # HDR excluded: diverged, often TE-rich, not evidence the outgroup has the sequence
os.makedirs(os.path.join(out_dir, "01_polarized"), exist_ok=True)


def find_pair(x, y):
    """-> (syri.out path, x_is_ref) for the pair in either orientation."""
    for d in (dir002, dir_cross):
        for a, b, x_ref in ((x, y, True), (y, x, False)):
            p = os.path.join(d, f"{a}_vs_{b}", f"{a}_vs_{b}.syri.out")
            if os.path.exists(p):
                return p, x_ref
    raise SystemExit(f"no syri output for {x} / {y}")


class Frame:
    """Alignment of frame genome F against outgroup O, indexed in F coordinates."""

    def __init__(self, f, o):
        path, f_ref = find_pair(f, o)
        ci = (0, 1, 2) if f_ref else (5, 6, 7)       # F chrom/start/end columns
        own_ins, o_ins = ("DEL", "INS") if f_ref else ("INS", "DEL")
        al = defaultdict(list)     # chrom -> (start, end, type)
        fonly = defaultdict(list)  # chrom -> (start, end): F-carried sequence absent in O
        oins = defaultdict(list)   # chrom -> (anchor, size): O-carried insertion
        for line in open(path):
            r = line.rstrip("\n").split("\t")
            t = r[10]
            if t in AL_TYPES and r[ci[0]] != "-":
                s, e = sorted((int(r[ci[1]]), int(r[ci[2]])))
                al[r[ci[0]]].append((s, e, t))
            elif t == own_ins:           # F has extra sequence: interval in F coords
                s, e = int(r[ci[1]]), int(r[ci[2]])
                if e - s >= SMALL_INDEL:
                    fonly[r[ci[0]]].append((s, e))
            elif t == o_ins:             # O has extra sequence: anchor in F coords
                osz = (int(r[7]) - int(r[6])) if f_ref else (int(r[2]) - int(r[1]))
                if osz >= SMALL_INDEL:
                    oins[r[ci[0]]].append((int(r[ci[1]]), osz))
        self.al = {c: self._index(sorted(v)) for c, v in al.items()}
        self.fonly = {c: self._merge(v) for c, v in fonly.items()}
        self.oins = {c: sorted(v) for c, v in oins.items()}

    @staticmethod
    def _index(blocks):
        pm, m = [], 0
        for b in blocks:
            m = max(m, b[1])
            pm.append(m)
        return blocks, [b[0] for b in blocks], pm

    @staticmethod
    def _merge(iv):
        out = []
        for s, e in sorted(iv):
            if out and s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        return out

    def _overlaps(self, chrom, s, e):
        if chrom not in self.al:
            return []
        blocks, starts, pm = self.al[chrom]
        i, res = bisect_right(starts, e) - 1, []
        while i >= 0 and pm[i] >= s:
            if blocks[i][1] >= s:
                res.append(blocks[i])
            i -= 1
        return res

    @staticmethod
    def _cover(ivs, s, e):
        """bp of [s,e] covered by the union of intervals ivs."""
        tot, cur = 0, s
        for a, b, *_ in sorted(ivs):
            a, b = max(a, cur), min(b, e)
            if b > a:
                tot += b - a
                cur = b
        return tot

    def has_fraction(self, chrom, s, e):
        """(aligned-and-not-F-only fraction, F-only-indel fraction) of F interval."""
        n = max(e - s, 1)
        cov = self._cover(self._overlaps(chrom, s, e), s, e)
        lost = self._cover([(a, b) for a, b in self.fonly.get(chrom, []) if b >= s and a <= e], s, e)
        return max(cov - lost, 0) / n, lost / n

    def call_carrier_F(self, chrom, s, e):
        has, lacks = self.has_fraction(chrom, s, e)
        if lacks >= FRAC:
            return "ins_F"       # O lacks the sequence F carries
        if has >= FRAC:
            return "del_G"       # O has it too -> G lost it
        return "amb"

    def call_carrier_G(self, chrom, anchor, size):
        near = [sz for a, sz in self.oins.get(chrom, []) if abs(a - anchor) <= ANCHOR_TOL]
        if any(size / SIZE_RATIO <= sz <= size * SIZE_RATIO for sz in near):
            return "del_F"       # O carries an insertion like G's -> F lost it
        if near:
            return "amb"
        has, lacks = self.has_fraction(chrom, anchor - FLANK, anchor + FLANK)
        if has >= FRAC and lacks == 0:
            return "ins_G"
        return "amb"

    def call_sv(self, chrom, s, e, t):
        bt = Counter()
        for a, b, bty in self._overlaps(chrom, s, e):
            bt[bty] += max(min(b, e) - max(a, s), 0)
        n = max(e - s, 1)
        syn = bt["SYNAL"] / n
        other = (bt["INVAL"] if t == "INV" else bt["TRANSAL"] + bt["INVTRAL"]) / n
        if syn >= SV_FRAC and syn > other:
            return "G"
        if other >= SV_FRAC and other > syn:
            return "F"
        return "amb"


def outgroups(a, b):
    ca, cb = a.startswith("Cbp"), b.startswith("Cbp")
    if ca and cb:
        return DIP[CBP_SET(a)], "cbp_cbp"
    s = CBP_SET(a) if ca else CBP_SET(b) if cb else SET_OF[a]
    return DIP["Cr" if s == "Co" else "Co"], ("dip_cbp" if ca or cb else "dip_dip")


def consensus(calls):
    c = Counter(x for x in calls if x != "amb")
    if not c:
        return "amb", 0, 0
    top = c.most_common()
    if len(top) > 1 and top[0][1] == top[1][1]:
        return "amb", sum(c.values()), 0
    return top[0][0], sum(c.values()), top[0][1]


frames = {}
def frame(f, o):
    if (f, o) not in frames:
        frames[(f, o)] = Frame(f, o)
    return frames[(f, o)]


def size_class(n):
    return ">=1kb" if n >= 1000 else "50bp-1kb"


# frame label -> lineage label: frame X's own lineage is "F", the other genome's is "G"
def to_lineage(call, x, other):
    if call == "amb":
        return call
    return call.replace("_F", "_" + x).replace("_G", "_" + other) if "_" in call else \
        (x if call == "F" else other)


def reconcile(cA, cB):
    """-> (final, confidence) from the A-frame and B-frame consensus."""
    if cA == cB:
        return cA, ("both" if cA != "amb" else "none")
    if cA == "amb":
        return cB, "B_frame_only"
    if cB == "amb":
        return cA, "A_frame_only"
    return "conflict", "conflict"


summ = Counter()
for line in open(pairs_tsv):
    if not line.strip():
        continue
    _, a, b = line.rstrip("\n").split("\t")
    outs, cls = outgroups(a, b)
    outs = [o for o in outs if o not in (a, b)]
    pair = f"{a}_vs_{b}"
    path = os.path.join(dir002, pair, pair + ".syri.out")
    fr = {a: [frame(a, o) for o in outs], b: [frame(b, o) for o in outs]}
    rows = []
    for line2 in open(path):
        r = line2.rstrip("\n").split("\t")
        t = r[10]
        ev = None
        if t == "DEL":          # A carries [r1,r2]; B has an anchor
            sz = int(r[2]) - int(r[1])
            ev = ("indel", "A", sz, (r[0], int(r[1]), int(r[2])), (r[5], int(r[6]), int(r[6])))
        elif t == "INS":        # B carries [r6,r7]; A has an anchor
            sz = int(r[7]) - int(r[6])
            ev = ("indel", "B", sz, (r[0], int(r[1]), int(r[1])), (r[5], int(r[6]), int(r[7])))
        elif t in ("INV", "TRANS", "INVTR", "DUP", "INVDP") and r[9] == "-":
            ia = (r[0], *sorted((int(r[1]), int(r[2]))))
            ib = (r[5], *sorted((int(r[6]), int(r[7]))))
            ev = (t, None, ia[2] - ia[1], ia, ib)
        if ev is None or ev[2] < MIN_SIZE:
            continue
        kind, carrier, size, ia, ib = ev
        if kind in ("DUP", "INVDP"):
            rows.append([pair, cls, kind, ia[0], ia[1], ia[2], size, "", "", "unpolarized", "", "", "", ""])
            summ[(pair, cls, a, b, kind, size_class(size), "unpolarized", "none")] += 1
            summ[(pair, cls, a, b, kind, size_class(size), "unpolarized", "none", "bp")] += size
            continue
        per = {}
        for x, other, mine, theirs in ((a, b, ia, ib), (b, a, ib, ia)):
            if kind == "indel":
                if carrier == ("A" if x == a else "B"):
                    cl = [f.call_carrier_F(*mine) for f in fr[x]]
                else:
                    cl = [f.call_carrier_G(mine[0], mine[1], size) for f in fr[x]]
            else:
                cl = [f.call_sv(*mine, kind) for f in fr[x]]
            per[x] = [to_lineage(c, "A" if x == a else "B", "B" if x == a else "A") for c in cl]
        cA, nA, _ = consensus(per[a])
        cB, nB, _ = consensus(per[b])
        final, conf = reconcile(cA, cB)
        rows.append([pair, cls, kind, ia[0], ia[1], ia[2], size, ";".join(per[a]),
                     ";".join(per[b]), final, conf, cA, cB, nA + nB])
        summ[(pair, cls, a, b, kind, size_class(size), final, conf)] += 1
        summ[(pair, cls, a, b, kind, size_class(size), final, conf, "bp")] += size
    with open(os.path.join(out_dir, "01_polarized", pair + ".tsv"), "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pair", "class", "kind", "chrA", "startA", "endA", "size",
                    "callsAframe(" + ",".join(outs) + ")", "callsBframe", "final", "confidence",
                    "consA", "consB", "n_calls"])
        w.writerows(rows)
    print(f"{pair}: {len(rows)} events, outgroups {outs}", flush=True)

keys = sorted({k[:8] for k in summ})
with open(os.path.join(out_dir, "polarization_summary.tsv"), "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["pair", "class", "A", "B", "kind", "size_class", "final", "confidence", "n", "bp"])
    for k in keys:
        w.writerow(list(k) + [summ[k], summ[k + ("bp",)]])
print("done", out_dir)
