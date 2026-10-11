import filecmp, os, sys
a, b = sys.argv[1], sys.argv[2]
skip = ('src_ts', 'run.log', 'plot.log', '__pycache__', 'compTime', 'profile.rank', 'bGlobal.txt')
fa = {f for f in os.listdir(a) if not f.startswith(skip) and os.path.isfile(os.path.join(a, f))}
fb = {f for f in os.listdir(b) if not f.startswith(skip) and os.path.isfile(os.path.join(b, f))}
ident = dateonly = 0; diff = []
for f in sorted(fa & fb):
    pa, pb = os.path.join(a, f), os.path.join(b, f)
    if filecmp.cmp(pa, pb, shallow=False): ident += 1; continue
    la, lb = open(pa, 'rb').read().split(b'\n'), open(pb, 'rb').read().split(b'\n')
    if len(la) == len(lb) and all(x == y or (b'# date =' in x and b'# date =' in y) for x, y in zip(la, lb)): dateonly += 1
    else: diff.append(f)
print('%s vs %s: %d files byte-identical, %d identical except the "# date =" header line, only-a %s, only-b %s, DIFFERING %s' % (os.path.basename(a), os.path.basename(b), ident, dateonly, sorted(fa - fb), sorted(fb - fa), diff))
