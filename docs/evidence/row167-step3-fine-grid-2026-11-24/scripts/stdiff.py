import sys,glob,os,numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R=sys.argv[1]; c=sys.argv[2]; out=sys.argv[3]
import re
def fl(x): return float(re.sub(r'(\d)([-+]\d{3})$',r'\1E\2',x))
def ld(p): return np.array([[fl(x) for x in l.split()] for l in open(p) if l.strip() and not l.lstrip().startswith('#') and l.split()[0][0] in '-0123456789.'])
names=sorted(os.path.basename(f)[:-4] for f in glob.glob(f'{R}/pre_{c}_dx250/*st*.txt') if 'bStations' not in f and 'correct' not in f)
rows=[]
for n in names:
    a=ld(f'{R}/pre_{c}_dx250/{n}.txt'); b=ld(f'{R}/post_{c}_dx250/{n}.txt')
    if a.shape!=b.shape: rows.append((n,'SHAPE',a.shape,b.shape)); continue
    S=max(np.abs(a[:,1]).max(),np.abs(b[:,1]).max(),1e-6)
    sr=max(np.abs(a[:,2]).max(),np.abs(b[:,2]).max(),1e-6)
    rows.append((n,np.abs(a[:,1]-b[:,1]).max()/S,np.abs(a[:,2]-b[:,2]).max()/sr,S,sr))
for r in rows: print(c,*[('%.4g'%x if isinstance(x,float) else x) for x in r])
ok=[r for r in rows if r[1]!='SHAPE']
print(c,'n=%d maxrel_slip=%.4g maxrel_sliprate=%.4g nonzero_slip_stations=%d'%(len(ok),max(r[1] for r in ok),max(r[2] for r in ok),sum(1 for r in ok if r[3]>1e-3)))
pick=[n for n in ('faultst000dp120','faultst040dp060','faultst080dp030') if n in names]
fig,ax=plt.subplots(len(pick),2,figsize=(9,2.2*len(pick)),squeeze=False)
for i,n in enumerate(pick):
    a=ld(f'{R}/pre_{c}_dx250/{n}.txt'); b=ld(f'{R}/post_{c}_dx250/{n}.txt')
    for j,(k,l) in enumerate(((1,'h-slip (m)'),(2,'h-slip rate (m/s)'))):
        ax[i,j].plot(a[:,0],a[:,k],'-',c='#0072B2',label='pre-fix'); ax[i,j].plot(b[:,0],b[:,k],'--',c='#D55E00',label='post-fix'); ax[i,j].set_title(n+' '+l,fontsize=8)
ax[0,0].legend(fontsize=7); fig.suptitle(c+' dx 250 m, same-name stations, no interpolation; cross-version consistency, not validation',fontsize=8); fig.tight_layout(); fig.savefig(out,dpi=110)
