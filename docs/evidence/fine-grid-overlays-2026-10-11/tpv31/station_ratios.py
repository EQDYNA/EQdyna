import numpy as np, os, glob, json, sys
H=os.path.expanduser('~/shared_dataset')
B={31:H+'/scec_cvws.tpv3132/raw/tpv31/barall.2-faultmod-25m-2014',32:H+'/scec_cvws.tpv3132/raw/tpv32/barall-faultmod-25m-2014',33:H+'/scec_cvws.tpv33/raw/tpv33/barall-faultmod-12.5m-2015'}
T={31:15,32:15,33:13}
def load(p):
    rows=[]
    for l in open(p):
        s=l.split()
        if not s or s[0].startswith('#'): continue
        try: rows.append([float(x.replace('E','e')) for x in s])
        except ValueError: continue
    return np.array(rows)
out={}
for c,b in B.items():
    r={'slip':[], 'psr':[], 'pv':[]}
    for f in sorted(os.listdir(b)):
        if f=='cplot': continue
        o=f'RUNS/tpv{c}_dx100/{f}'
        if not os.path.exists(o): o+='.txt'
        if not os.path.exists(o): continue
        a,z=load(o),load(b+'/'+f); a=a[a[:,0]<=T[c]+1e-6]; z=z[z[:,0]<=T[c]+1e-6]
        if f.startswith('fault'):
            sa,sz=np.hypot(a[-1,1],a[-1,4]),np.hypot(z[-1,1],z[-1,4])
            if sz>0.01: r['slip'].append(sa/sz)
            pa,pz=np.hypot(a[:,2],a[:,5]).max(),np.hypot(z[:,2],z[:,5]).max()
            if pz>0.01: r['psr'].append(pa/pz)
        else:
            pa,pz=np.abs(a[:,[2,4,6]]).max(),np.abs(z[:,[2,4,6]]).max()
            if pz>1e-3: r['pv'].append(pa/pz)
    out[c]={k:dict(n=len(v),median=float(np.median(v)),min=float(np.min(v)),max=float(np.max(v))) for k,v in r.items()}
    print(c,json.dumps({k:'n%d med %.3f [%.3f-%.3f]'%(d['n'],d['median'],d['min'],d['max']) for k,d in out[c].items()}))
    json.dump(out[c],open(f'docs/evidence/fine-grid-overlays-2026-10-11/tpv{c}/station_ratios_faultmod.json','w'),indent=1)
