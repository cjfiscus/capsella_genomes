#!/usr/bin/env python3
"""Collect SyRI results for all pairs in data/pairs.tsv into three tables.

usage: 003a_summarize.py <pairs.tsv> <02_syri dir> <01_norm dir> <out dir>

  sv_summary.tsv        every row of every <pair>.syri.summary (long format)
  pair_metrics.tsv      one row per pair: headline numbers + fraction of each genome
                        that is syntenic / not aligned
  sv_size_filtered.tsv  counts and bp of SVs >= 50 bp and >= 1 kb, straight from
                        syri.out (blunts the preset sensitivity of raw counts,
                        which are inflated by alignment fragmentation)

Fails loudly if any pair is missing or empty.
"""
import csv
import os
import sys

pairs_tsv, syri_dir, norm_dir, out_dir = sys.argv[1:5]
MIN_SIZES = (50, 1000)
os.makedirs(out_dir, exist_ok=True)

pairs = [l.rstrip("\n").split("\t") for l in open(pairs_tsv) if l.strip()]


def genome_len(name):
    fai = os.path.join(norm_dir, name + ".fa.fai")
    with open(fai) as fh:
        return sum(int(l.split("\t")[1]) for l in fh)


def num(x):
    return None if x in ("-", "") else int(x)


# ---- 1) long summary + per-pair metrics --------------------------------------
KEYS = {
    "Syntenic regions": "syn",
    "Inversions": "inv",
    "Translocations": "trans",
    "Duplications (reference)": "dup_ref",
    "Duplications (query)": "dup_qry",
    "Not aligned (reference)": "notal_ref",
    "Not aligned (query)": "notal_qry",
    "SNPs": "snp",
    "Insertions": "ins",
    "Deletions": "del",
    "Highly diverged": "hdr",
}

long_rows, metric_rows = [], []
for set_, ref, qry in pairs:
    pair = f"{ref}_vs_{qry}"
    path = os.path.join(syri_dir, pair, f"{pair}.syri.summary")
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        sys.exit(f"ERROR: missing/empty {path}")
    block, parsed = None, {}
    for line in open(path):
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith("#Structural"):
            block = "structural"
        elif line.startswith("#Sequence"):
            block = "sequence"
        elif line.startswith("#"):
            continue
        else:
            typ, count, lref, lqry = line.split("\t")
            long_rows.append([pair, set_, ref, qry, block, typ, int(count), num(lref), num(lqry)])
            if typ in KEYS:
                parsed[KEYS[typ]] = (int(count), num(lref), num(lqry))
    missing = [k for k in KEYS.values() if k not in parsed]
    if missing:
        sys.exit(f"ERROR: {path} lacks expected rows: {missing}")
    rlen, qlen = genome_len(ref), genome_len(qry)
    metric_rows.append([
        pair, set_, ref, qry, rlen, qlen,
        parsed["syn"][1], round(parsed["syn"][1] / rlen, 4), parsed["syn"][2], round(parsed["syn"][2] / qlen, 4),
        parsed["notal_ref"][1], round(parsed["notal_ref"][1] / rlen, 4),
        parsed["notal_qry"][2], round(parsed["notal_qry"][2] / qlen, 4),
        parsed["inv"][0], parsed["trans"][0], parsed["dup_ref"][0], parsed["dup_qry"][0],
        parsed["snp"][0], parsed["ins"][0], parsed["del"][0],
    ])

with open(os.path.join(out_dir, "sv_summary.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["pair", "set", "ref", "qry", "block", "type", "count", "len_ref", "len_qry"])
    w.writerows([["NA" if v is None else v for v in r] for r in long_rows])

with open(os.path.join(out_dir, "pair_metrics.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["pair", "set", "ref", "qry", "ref_len", "qry_len",
                "syn_ref_bp", "syn_ref_frac", "syn_qry_bp", "syn_qry_frac",
                "notal_ref_bp", "notal_ref_frac", "notal_qry_bp", "notal_qry_frac",
                "n_inv", "n_trans", "n_dup_ref", "n_dup_qry", "n_snp", "n_ins", "n_del"])
    w.writerows(metric_rows)

# ---- 2) size-filtered SV counts from syri.out --------------------------------
# syri.out cols: 1 refchr 2 refstart 3 refend 4 refseq 5 qryseq 6 qrychr 7 qrystart 8 qryend
#                9 ID 10 parent 11 type 12 copystatus
# structural (top-level, parent '-'): INV TRANS INVTR DUP INVDP NOTAL
# small-variant children: INS (size = len(qryseq)-1), DEL (size = len(refseq)-1); the -1 is the anchor base
STRUCT = ("INV", "TRANS", "INVTR", "DUP", "INVDP")
TYPES = STRUCT + ("NOTAL_ref", "NOTAL_qry", "INS", "DEL")

size_rows = []
for set_, ref, qry in pairs:
    pair = f"{ref}_vs_{qry}"
    path = os.path.join(syri_dir, pair, f"{pair}.syri.out")
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        sys.exit(f"ERROR: missing/empty {path}")
    n = {(t, m): 0 for t in TYPES for m in MIN_SIZES}
    bp = {(t, m): 0 for t in TYPES for m in MIN_SIZES}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            typ = f[10]
            if typ in ("INS", "DEL"):
                size = (len(f[4]) if typ == "INS" else len(f[3])) - 1
                key = typ
            elif f[9] == "-" and typ in STRUCT:
                rs = abs(int(f[2]) - int(f[1])) + 1
                qs = abs(int(f[7]) - int(f[6])) + 1
                size, key = max(rs, qs), typ
            elif f[9] == "-" and typ == "NOTAL":
                if f[0] != "-":
                    size, key = abs(int(f[2]) - int(f[1])) + 1, "NOTAL_ref"
                else:
                    size, key = abs(int(f[7]) - int(f[6])) + 1, "NOTAL_qry"
            else:
                continue
            for m in MIN_SIZES:
                if size >= m:
                    n[(key, m)] += 1
                    bp[(key, m)] += size
    for t in TYPES:
        for m in MIN_SIZES:
            size_rows.append([pair, set_, ref, qry, t, m, n[(t, m)], bp[(t, m)]])
    print(f"[{pair}] size-filtered counts done", flush=True)

with open(os.path.join(out_dir, "sv_size_filtered.tsv"), "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["pair", "set", "ref", "qry", "type", "min_size", "count", "bp"])
    w.writerows(size_rows)

print(f"wrote {len(long_rows)} summary rows, {len(metric_rows)} pairs, {len(size_rows)} size rows -> {out_dir}")
