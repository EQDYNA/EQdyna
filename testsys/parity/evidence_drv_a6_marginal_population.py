"""Item 32: isolate the PORT's contribution to drv.a6 rupture-arrival flips.

Against the committed 4-rank reference, the serial python backends show
415-423 flips. Fortran run SERIALLY against that same 4-rank reference shows
329 -- with Fortran on both sides -- so most of what the port is charged with
is decomposition, not port difference.

This compares serial python against a SERIAL FORTRAN reference: same
decomposition on both sides, so the decomposition term is zero by construction
and what remains is the port alone.
"""
import sys, numpy as np
sys.path.insert(0, '/home/utig5/dliu/EQdyna')
from testsys import compare, matrix, frt_canonical

S = '/tmp/claude-16759/-home-utig5-dliu-EQdyna/e618789a-38a5-4235-9ca3-c8a3a0fe9ba6/scratchpad'
REF_SERIAL = S + '/drv_serial_fortran_ref.frt'
CASE = 'test.drv.a6'

def flips(a_path, b_path, label):
    a, b = compare.align_two_frt_files(a_path, b_path)
    d = compare.flip_decomposition(a, b)
    print('  %-34s flips %4d/%d  (ref-only %3d, run-only %3d, timing-shift %3d)  '
          'median|dfnft| %.4f s  phys_max %.3e'
          % (label, d['total_flips'], matrix.DRV_A6['total_flip_bound'],
             d['n_only_ref'], d['n_only_run'], d['n_timing_shifts'],
             d['median_fnft_diff'], d['phys_max_diff']))
    return d['total_flips']

print('drv.a6 rupture-arrival flips, 5151 fault nodes\n')
print('A. against the COMMITTED 4-RANK reference (what the gate uses):')
ref4 = compare.reference_path(CASE)
f4 = {}
for b in ('numpy', 'jax'):
    p = '%s/drv_py_%s/test.drv.a6/frt.txt0' % (S, b)
    try: f4[b] = flips(ref4, p, 'python-%s vs 4-rank ref' % b)
    except Exception as e: print('  python-%s: %s' % (b, e))
try: f4['fortran'] = flips(ref4, REF_SERIAL, 'fortran-SERIAL vs 4-rank ref')
except Exception as e: print('  fortran-serial: %s' % e)

print('\nB. against the SERIAL FORTRAN reference (decomposition term removed):')
fs = {}
for b in ('numpy', 'jax'):
    p = '%s/drv_py_%s/test.drv.a6/frt.txt0' % (S, b)
    try: fs[b] = flips(REF_SERIAL, p, 'python-%s vs fortran-serial' % b)
    except Exception as e: print('  python-%s: %s' % (b, e))

print('\nVERDICT')
if 'fortran' in f4:
    print('  decomposition alone (fortran both sides) : %4d flips' % f4['fortran'])
for b in ('numpy', 'jax'):
    if b in f4 and b in fs:
        print('  python-%-5s  vs 4-rank %4d  ->  vs serial %4d   '
              'port-only share %.0f%%' % (b, f4[b], fs[b], 100.0*fs[b]/max(f4[b],1)))
