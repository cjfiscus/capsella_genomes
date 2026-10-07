#!/usr/bin/env python3
"""Diploid-diploid baseline divergence (ANALYSES_INVENTORY B1 #9), from 002_alignment output.
Per pair: SNP/kb, indel/kb (syri INS+DEL rows within alignments), and large-SV/Mb
over syntenic reference bp; classes = diploid-diploid, Cbp-vs-diploid, Cbp-vs-Cbp.
usage: python3 scripts/000_baseline.py   (run from this analysis dir)"""
import csv, statistics as st, collections
SRC = '../002_alignment/results/03_summary/pair_metrics.tsv'
dip = lambda g: not g.startswith('Cbp')
rows = []
for r in csv.DictReader(open(SRC), delimiter='\t'):
    kb = float(r['syn_ref_bp']) / 1e3
    cls = ('diploid-diploid' if dip(r['ref']) and dip(r['qry']) else
           'Cbp-Cbp' if not dip(r['ref']) and not dip(r['qry']) else 'Cbp-diploid')
    sv = sum(int(r[k]) for k in ('n_inv', 'n_trans', 'n_dup_ref'))
    rows.append(dict(pair=r['pair'], set=r['set'], cls=cls,
        snp_per_kb=int(r['n_snp'])/kb, indel_per_kb=(int(r['n_ins'])+int(r['n_del']))/kb,
        sv_per_mb=sv/(kb/1e3), syn_frac=min(float(r['syn_ref_frac']), float(r['syn_qry_frac']))))
f = ['pair', 'set', 'cls', 'snp_per_kb', 'indel_per_kb', 'sv_per_mb', 'syn_frac']
with open('results/pair_baseline.tsv', 'w') as o:
    w = csv.DictWriter(o, f, delimiter='\t'); w.writeheader()
    for r in rows: w.writerow({k: (f'{v:.3f}' if isinstance(v, float) else v) for k, v in r.items()})
g = collections.defaultdict(list)
for r in rows: g[(r['set'], r['cls'])].append(r)
with open('results/class_summary.tsv', 'w') as o:
    o.write('set\tclass\tn\t' + '\t'.join(f'{m}_med\t{m}_min\t{m}_max' for m in f[3:]) + '\n')
    for (s, c), v in sorted(g.items()):
        o.write(f'{s}\t{c}\t{len(v)}')
        for m in f[3:]:
            x = [i[m] for i in v]; o.write(f'\t{st.median(x):.3f}\t{min(x):.3f}\t{max(x):.3f}')
        o.write('\n')
