#!/usr/bin/env python3
"""TE fraction of syri regions for every pair whose two genomes both have a
results/04_te/<norm>.te.bed (made by 004a_prep_te.py; coordinates = 01_norm).
Writes results/04_te/te_overlap.tsv: pair, type, side, Mb, TE_Mb, TE_frac.
Genomes without a TE bed (dropped from 001_edta) are skipped."""
import bisect, collections, os
D = os.path.dirname(os.path.abspath(__file__)) + '/../results'
def load(p):
    d = collections.defaultdict(list)
    for l in open(p):
        c, s, e = l.split()[:3]; d[c].append((int(s) + 1, int(e)))
    return d
def cov(m, c, s, e):
    iv = m.get(c, []); i = max(0, bisect.bisect_left(iv, (s, 0)) - 1); tot = 0
    while i < len(iv) and iv[i][0] <= e:
        tot += max(0, min(e, iv[i][1]) - max(s, iv[i][0]) + 1); i += 1
    return tot
TYPES = ('NOTAL', 'HDR', 'SYN', 'INV', 'DUP', 'TRANS', 'INVTR', 'INVDP')
cache = {}
def te(n):
    if n not in cache: cache[n] = load(f'{D}/04_te/{n}.te.bed')
    return cache[n]
with open(f'{D}/04_te/te_overlap.tsv', 'w') as o:
    o.write('pair\ttype\tside\tMb\tTE_Mb\tTE_frac\n')
    for l in open(f'{D}/../data/pairs.tsv'):
        f = l.split()
        r, q = f[1], f[2]
        if not all(os.path.exists(f'{D}/04_te/{x}.te.bed') for x in (r, q)): continue
        pair = f'{r}_vs_{q}'; M = {'ref': te(r), 'qry': te(q)}
        acc = collections.defaultdict(lambda: [0, 0])
        for ln in open(f'{D}/02_syri/{pair}/{pair}.syri.out'):
            x = ln.rstrip('\n').split('\t')
            if x[10] not in TYPES: continue
            for side, c, s, e in (('ref', *x[0:3]), ('qry', *x[5:8])):
                if c == '-': continue
                s, e = sorted((int(s), int(e)))
                a = acc[(x[10], side)]; a[0] += e - s + 1; a[1] += cov(M[side], c, s, e)
        for (t, side), (b, n) in sorted(acc.items()):
            o.write(f'{pair}\t{t}\t{side}\t{b/1e6:.2f}\t{n/1e6:.2f}\t{n/b:.3f}\n')
