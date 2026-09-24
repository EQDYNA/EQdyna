#! /usr/bin/env python3
"""Regression guard: box tenancy is measurable on a host with no NUMA
topology, and a profile-capture SystemExit fails one cell, not the sweep.

Incident (2026-09-23, PR #6 CI smoke): the run-profile record calls
ledger.box_tenancy on every cell; on the GitHub runner `numactl --hardware`
reports no topology, box_tenancy raised SystemExit, and run_e2e's
`except Exception` did not catch it -- the smoke process exited 1 after the
first cell. Two checks, each mutation-tested:
  1. with numa_topology() stubbed to {} (the runner's answer), box_tenancy
     returns total == the number of cpuN lines in /proc/stat;
  2. _run_profile_capture turns a SystemExit from capture_run into
     (False, ['PROFILE CAPTURE FAILED ...']) instead of propagating it.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'testsys'))
sys.path.insert(0, os.path.join(ROOT, 'testsys', 'perf'))
sys.path.insert(0, os.path.join(ROOT, 'testsys', 'e2e'))

import ledger  # noqa: E402
import run_numa_scaling  # noqa: E402
import run_e2e  # noqa: E402


def main():
    fails = []
    n_proc = sum(1 for ln in open('/proc/stat')
                 if ln.startswith('cpu') and ln[3:4].isdigit())
    real = run_numa_scaling.numa_topology
    run_numa_scaling.numa_topology = lambda *a, **k: {}
    try:
        t = ledger.box_tenancy(ledger.TENANCY_REFERENCE_CEILING)
        if t.get('total') != n_proc:
            fails.append('no-NUMA box_tenancy total=%r, /proc/stat lists %d cpus'
                         % (t.get('total'), n_proc))
    except SystemExit as exc:
        fails.append('no-NUMA box_tenancy raised SystemExit: %s' % exc)
    finally:
        run_numa_scaling.numa_topology = real

    real_cap = run_e2e.profile_record.capture_run

    def boom(*a, **k):
        raise SystemExit('FAIL: simulated tenancy failure')
    run_e2e.profile_record.capture_run = boom
    try:
        ok, lines = run_e2e._run_profile_capture('/nonexistent', 'test.tpv8',
                                                 'fortran', '5.0', 'deadbee', False)
        if ok or not lines or not lines[0].startswith('PROFILE CAPTURE FAILED'):
            fails.append('_run_profile_capture on SystemExit returned %r' % ((ok, lines),))
    except SystemExit as exc:
        fails.append('_run_profile_capture let SystemExit escape: %s' % exc)
    except TypeError as exc:
        fails.append('_run_profile_capture signature changed, update this guard: %s' % exc)
    finally:
        run_e2e.profile_record.capture_run = real_cap

    if fails:
        print('FAIL test_tenancy_without_numa:\n  - ' + '\n  - '.join(fails))
        return 1
    print('SUCCESS test_tenancy_without_numa (no-NUMA tenancy total=%d; '
          'SystemExit contained to one cell)' % n_proc)
    return 0


if __name__ == '__main__':
    sys.exit(main())
