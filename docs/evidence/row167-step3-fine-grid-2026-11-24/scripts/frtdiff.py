import sys,glob,numpy as np
def ld(d):
    rows=[]
    for f in glob.glob(d+'/frt.txt*'):
        for l in open(f):
            t=l.split()
            if len(t)>=13:
                try: rows.append([float(x) for x in t])
                except: pass
    a=np.array(rows); return a[np.lexsort((a[:,2],a[:,1],a[:,0]))] if False else a
c,dx=sys.argv[1:3]
a=ld(f'runs/pre_{c}_dx{dx}'); b=ld(f'runs/post_{c}_dx{dx}')
k=lambda m:np.lexsort((np.round(m[:,2],1),np.round(m[:,1],1),np.round(m[:,0],1)))
a=a[k(a)];b=b[k(b)]
print(c,dx,a.shape,b.shape)
if a.shape==b.shape:
    d=np.abs(a-b); print('cols with diff:',[ (i+1,float(d[:,i].max())) for i in range(d.shape[1]) if d[:,i].max()>0])
    r=(d[:,10:13].max(1)>0); print('rows changed (cols11-13):',int(r.sum()),'max depth-coord of changed rows',float(np.abs(a[r,2]).max()) if r.any() else None)
