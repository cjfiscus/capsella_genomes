import csv, glob, os, collections
MIN=1_000_000
pairs=[l.split() for l in open('data/pairs.tsv')]
# per genome: list of (chr,start,end,partner)
G=collections.defaultdict(list)
for s,r,q in pairs:
    f=glob.glob(f'results/02_syri/{r}_vs_{q}/*.syri.out')[0]
    for l in open(f):
        c=l.rstrip('\n').split('\t')
        if c[10]!='INV': continue
        rs,re_=int(c[1]),int(c[2]); 
        if re_-rs+1>=MIN: G[r].append((c[0],rs,re_,q))
        qs,qe=sorted((int(c[6]),int(c[7])))
        if qe-qs+1>=MIN: G[q].append((c[5],qs,qe,r))
sets={s:set() for s,_,_ in pairs}
for s,r,q in pairs: sets[s]|={r,q}
gset={g:s for s,v in sets.items() for g in v}
print("genome\tset\tchr\tstart\tend\tlen_Mb\tn_partners_supporting/n_partners\tsupporters")
for g in sorted(G):
    recs=sorted(G[g]); clusters=[]
    for c,s,e,p in recs:
        for cl in clusters:
            if cl['c']==c and min(e,cl['e'])-max(s,cl['s'])>0.5*min(e-s,cl['e']-cl['s']):
                cl['s']=min(cl['s'],s);cl['e']=max(cl['e'],e);cl['p'].add(p);break
        else: clusters.append(dict(c=c,s=s,e=e,p={p}))
    n=len(sets[gset[g]])-1
    for cl in clusters:
        print(f"{g}\t{gset[g]}\t{cl['c']}\t{cl['s']}\t{cl['e']}\t{(cl['e']-cl['s'])/1e6:.1f}\t{len(cl['p'])}/{n}\t{','.join(sorted(cl['p']))}")
