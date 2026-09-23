"""Aggregate the three A/B repetitions. Paired WITHIN a repetition: each arm
is divided by the base arm measured in the SAME repetition on the SAME cpu
set, so box drift between repetitions cancels instead of being averaged in."""
import json, glob, statistics as st

files = sorted(glob.glob('docs/perf_snapshots/jaxmpi_ab_*.json'))
print('repetitions:', len(files))
for f in files:
    print('  ', f)
data = {}
sets = {}
for i, f in enumerate(files):
    m = json.load(open(f))
    for row in m['rows']:
        n = row['ranks']
        sets[(i, n)] = row['cpus']
        for pos, (arm, r) in enumerate(row['arms'].items()):
            if r.get('failed'):
                continue
            data.setdefault((n, arm), []).append(
                (i, r['ms_per_step'], r['min_eff'], pos, r['mpi_ms'],
                 r['wait_ms'], r['rank_ms']))

print('\nRAW ms/step by repetition (arm ORDER was shuffled between them)')
print('%-6s %-5s %9s %9s %9s  %9s %7s' %
      ('ranks', 'arm', 'rep1', 'rep2', 'rep3', 'mean', 'sd%'))
for n in (4, 8, 12):
    for arm in ('base', 'b', 'a', 'ab'):
        v = sorted(data[(n, arm)])
        ms = [x[1] for x in v]
        print('%-6d %-5s %9.2f %9.2f %9.2f  %9.2f %6.1f%%'
              % (n, arm, ms[0], ms[1], ms[2], st.mean(ms),
                 100 * st.pstdev(ms) / st.mean(ms)))

print('\nPAIRED speedup vs the base arm of the SAME repetition')
print('%-6s %-5s %8s %8s %8s  %8s %8s' %
      ('ranks', 'arm', 'rep1', 'rep2', 'rep3', 'mean', 'sd'))
sp = {}
for n in (4, 8, 12):
    base = {x[0]: x[1] for x in data[(n, 'base')]}
    for arm in ('b', 'a', 'ab'):
        r = [base[p] / ms for p, ms, *_ in sorted(data[(n, arm)])]
        sp[(n, arm)] = r
        print('%-6d %-5s %8.3f %8.3f %8.3f  %8.3f %8.3f'
              % (n, arm, r[0], r[1], r[2], st.mean(r), st.pstdev(r)))

print('\nINTERACTION ab / (a * b) -- 1.0 means the two changes are additive')
for n in (4, 8, 12):
    v = [sp[(n, 'ab')][i] / (sp[(n, 'a')][i] * sp[(n, 'b')][i]) for i in range(3)]
    print('  %2d ranks  %s  mean %.3f'
          % (n, ' '.join('%.3f' % x for x in v), st.mean(v)))

print('\nMEASURED NOISE FLOOR: base-arm spread across repetitions on the same')
print('rank count (same case, same code, same cpu policy)')
for n in (4, 8, 12):
    ms = [x[1] for x in data[(n, 'base')]]
    print('  %2d ranks  %s  spread %.1f%% of mean'
          % (n, ' '.join('%.2f' % x for x in ms),
             100 * (max(ms) - min(ms)) / st.mean(ms)))

print('\ncpu set per repetition / rank count')
for k in sorted(sets):
    print('  rep%d %2dr %s' % (k[0] + 1, k[1], sets[k]))

print('\nexchange ms/step (mean over ranks) | barrier wait ms/step (mean)')
for n in (4, 8, 12):
    for arm in ('base', 'b', 'a', 'ab'):
        v = sorted(data[(n, arm)])
        ex = [st.mean(x[4]) for x in v]
        wt = [st.mean(x[5]) for x in v]
        print('  %2dr %-5s ex %s | wait %s' % (
            n, arm, ' '.join('%6.2f' % x for x in ex),
            ' '.join('%6.2f' % x for x in wt)))

print('\nmin EFFECTIVE_CORES over every rank of every point')
allm = [x[2] for v in data.values() for x in v]
print('  %d points, min %.2f, all >= 0.95: %s'
      % (len(allm), min(allm), all(m >= 0.95 for m in allm)))

