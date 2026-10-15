#!/usr/bin/env python3
"""Build per-01_norm-genome TE interval files (results/04_te/<norm>.te.bed, merged,
0-based BED) from the 001_edta panlib TEanno.gff3, in 01_norm coordinates.
Steps: decode EDTA '_J...' seqids (in-order match of X.fa.mod to X.fa; lengths verified),
keep only chromosomes in the chrmap for that subgenome, flip '-' strand chromosomes
(chrmap strand) using the chromosome length. Cbp2-2 uses the whole-genome gff for both subgenomes."""
import csv, collections, os, sys
E = '/rhome/cfisc004/bigdata/capsella_genomes/analyses/001_edta/results'
R = os.path.dirname(os.path.abspath(__file__)) + '/../results'
# norm genome -> (gff, .fa, .fa.mod, chrmap genome, subgenome or None)
P = E + '/02_panlib/'
SRC = {  # chrmap name -> (dir prefix)
 'Cbp2-2': P + 'Cbp2-2', 'CbpRussian': P + 'CbpRussian', 'Co39': P + 'Co39', 'Cr145': P + 'Cr145',
 'CoMosher': P + 'CoMosher', 'CrMosher': P + 'CrMosher',
 'CgMosher': E + '/03_cgmosher_renamed/CgMosher'}
def lens(fa):
    d, n = [], None
    for l in open(fa):
        if l[0] == '>': d.append([l[1:].split()[0], 0])
        else: d[-1][1] += len(l.strip())
    return d
def merge(v):
    v.sort(); out = []
    for s, e in v:
        if out and s <= out[-1][1]: out[-1][1] = max(out[-1][1], e)
        else: out.append([s, e])
    return out
for g, pre in SRC.items():
    a, b = lens(pre + '.fa'), lens(pre + '.fa.mod')
    assert len(a) == len(b) and all(x[1] == y[1] for x, y in zip(a, b)), g
    dec = {y[0]: x[0] for x, y in zip(a, b)}; length = {x[0]: x[1] for x in a}
    cm = list(csv.DictReader(open(f'{R}/01_norm/{g}.chrmap.tsv'), delimiter='\t'))
    iv = collections.defaultdict(list)
    for l in open(pre + '.fa.mod.EDTA.TEanno.gff3'):
        if l[0] == '#': continue
        f = l.split('\t')
        if len(f) < 5: continue
        iv[dec.get(f[0], f[0])].append((int(f[3]) - 1, int(f[4])))
    for sub in sorted({r['subgenome'] for r in cm}) if g.startswith('Cbp') else [None]:
        name = g if sub is None else f'{g}{sub}'
        n = 0
        with open(f'{R}/04_te/{name}.te.bed', 'w') as o:
            for r in cm:
                if sub and r['subgenome'] != sub: continue
                c = r['qchrom']; L = length[c]
                m = merge(list(iv.get(c, [])))
                if r['strand'] == '-': m = sorted([L - e, L - s] for s, e in m)
                for s, e in m: o.write(f'{c}\t{s}\t{e}\n')
                n += sum(e - s for s, e in m)
        print(name, f'{n/1e6:.1f} Mb TE', file=sys.stderr)
