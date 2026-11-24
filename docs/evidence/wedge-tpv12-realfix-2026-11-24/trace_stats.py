import sys, numpy as np
rows=[]
for l in open(sys.argv[1]):
    p=l.split()
    try: rows.append([float(x) for x in p])
    except ValueError: pass
d=np.array(rows); t=d[:,0]; n=d[:,7]   # column #8 = normal stress (MPa), + = tension
i=n.argmax(); k=np.abs(t-4.0).argmin()
print(f"{sys.argv[1].split('/')[-2]}: samples {len(t)} t_end {t[-1]:.3f}  n_max {n[i]:+.4f} MPa at t={t[i]:.3f}  n(4s) {n[k]:+.4f}  n_min {n.min():+.4f}  samples>0.1MPa {(n>0.1).sum()}  >1MPa {(n>1).sum()}  last-slip(col6?) see file")
