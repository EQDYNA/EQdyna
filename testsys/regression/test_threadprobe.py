#! /usr/bin/env python3
"""
Regression guard for testsys/perf/threadprobe.py (pathway item 47a: item 47's
probe was never committed and was lost; this is its rewrite). Behaviour only,
no box state, no subprocess (rule 10a):
  - summarize() computes spread and straggler ratios from per-worker timings;
  - it RAISES on a non-positive timing and on an empty result set (rule 2),
    rather than reporting a number.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'perf'))
import threadprobe as tp  # noqa: E402


def main():
    fails = []
    rows = [dict(wall_s=1.0, cpu_s=1.0, iters=100), dict(wall_s=2.0, cpu_s=1.0, iters=100),
            dict(wall_s=3.0, cpu_s=2.97, iters=100)]
    s = tp.summarize(rows)
    if (round(s['ms_min'], 6), round(s['ms_max'], 6), round(s['spread_max_over_min'], 6),
            round(s['straggler_max_over_mean'], 6)) != (10.0, 30.0, 3.0, 1.5):
        fails.append('summarize ratios wrong: %r' % s)
    if [round(e, 2) for e in s['effective_cores']] != [1.0, 0.5, 0.99]:
        fails.append('EFFECTIVE_CORES wrong: %r' % s['effective_cores'])
    ov = tp.summarize([dict(wall_s=1.0, cpu_s=1.0, iters=10, t_start=0.0, t_end=1.0),
                       dict(wall_s=1.0, cpu_s=1.0, iters=10, t_start=0.5, t_end=1.5)])['overlap']
    if round(ov, 6) != 0.5:
        fails.append('overlap of two half-overlapping windows is %r, want 0.5' % ov)
    ov2 = tp.summarize([dict(wall_s=1.0, cpu_s=1.0, iters=10, t_start=0.0, t_end=1.0),
                        dict(wall_s=2.68, cpu_s=2.6, iters=10, t_start=0.0, t_end=2.68)])['overlap']
    if round(ov2, 6) != 1.0:
        fails.append('a real 2.68x straggler that started together reads overlap %r, '
                     'want 1.0 -- the gate would refuse the result it exists to measure' % ov2)
    for bad, label in (([], 'empty'), ([dict(wall_s=0.0, cpu_s=0.0, iters=10)], 'zero wall'),
                       ([dict(wall_s=-1.0, cpu_s=0.0, iters=10)], 'negative wall')):
        try:
            tp.summarize(bad)
            fails.append('summarize(%s) did not raise' % label)
        except ValueError:
            pass
        except Exception as e:  # any other error = the refusal is missing
            fails.append('summarize(%s) raised %r instead of refusing' % (label, e))
    if fails:
        print('FAIL test_threadprobe')
        for f in fails:
            print('  -', f)
        return 1
    print('SUCCESS test_threadprobe (spread 3.00x, straggler 1.50x on a 10/20/30 ms '
          'fixture; empty, zero and negative timings refused)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
