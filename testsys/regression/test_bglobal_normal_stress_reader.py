#! /usr/bin/env python3
"""
Regression guard (#189/#200): the REAL Fortran readglobal and the REAL Python
read_bglobal give the same verdict on bGlobal.txt's last two lines,
"C_fric_pc_L fric_pc_L" and "C_normal_stress_caps min_norm max_norm".

One tiny planar friclaw-4 case (test_rsf_noncompressive_stop.py's control,
1 rank, 3 steps) is set up once; then each variant below rewrites those two
lines and runs bin/eqdyna and `python3 -m eqdyna ... --backend numpy`:

  - the file case.setup writes                      -> both exit 0;
  - an invalid caps switch, or not max_norm <= min_norm < 0, or a set
    fric_pc_L that is NaN / +-inf                   -> both exit 10
                                                      (ERR_CFG_NORMAL_STRESS_INPUTS_INVALID);
  - a stale "set L" line (old one-number line, set flag 2), or a stale
    caps line (missing, short)                      -> Fortran exits 22
                                                      (ERR_INPUT_FILE_STALE,
                                                      readglobal's two
                                                      stopStaleGlobal calls);
                                                      python fails with the
                                                      same "older case.setup"
                                                      verdict.

Departure, stated (it predates these lines): read_bglobal's stale() raises a
plain ValueError for every stale bGlobal.txt line, so python exits 1 there,
not 22. This test pins the verdict text on the python side, the code on the
Fortran side.

test_nstress_sign_convention.py drives read_bglobal alone on the same
fixtures; this one is the Fortran side and the cross-code agreement.

Fails (never skips) when mpirun or bin/eqdyna is missing.
"""
import os
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_rsf_noncompressive_stop as T  # noqa: E402  (the tiny case + runners)

ERR_CFG_NORMAL_STRESS_INPUTS_INVALID = 10
ERR_INPUT_FILE_STALE = 22
STALE_TEXT = 'older case.setup'

# (label, last two lines of bGlobal.txt, expected Fortran exit code; python
# exits the same except on a stale file, see above)
VARIANTS = [
    ('as case.setup writes it', '0 0.0\n-1 -10000000.0 -40000000.0\n', 0),
    ('caps switch 2', '0 0.0\n2 -1e7 -4e7\n', ERR_CFG_NORMAL_STRESS_INPUTS_INVALID),
    ('caps min_norm < max_norm', '0 0.0\n-1 -4e7 -1e7\n', ERR_CFG_NORMAL_STRESS_INPUTS_INVALID),
    ('set fric_pc_L nan', '1 nan\n-1 -1e7 -4e7\n', ERR_CFG_NORMAL_STRESS_INPUTS_INVALID),
    ('set fric_pc_L inf', '1 inf\n-1 -1e7 -4e7\n', ERR_CFG_NORMAL_STRESS_INPUTS_INVALID),
    ('set fric_pc_L -inf', '1 -inf\n-1 -1e7 -4e7\n', ERR_CFG_NORMAL_STRESS_INPUTS_INVALID),
    ('stale "set L": old one-number line', '0.004\n-1 -1e7 -4e7\n', ERR_INPUT_FILE_STALE),
    ('stale "set L": set flag 2', '2 0.1\n-1 -1e7 -4e7\n', ERR_INPUT_FILE_STALE),
    ('stale caps line: missing', '0 0.0\n', ERR_INPUT_FILE_STALE),
    ('stale caps line: short', '0 0.0\n-1 -1e7\n', ERR_INPUT_FILE_STALE),
]


def _base_text(path):
    """bGlobal.txt without its last two lines."""
    lines = open(path).read().splitlines(keepends=True)
    tail = ''.join(lines[-2:])
    if tail != '0 0.0\n-1 -10000000.0 -40000000.0\n':
        raise AssertionError('case.setup wrote an unexpected bGlobal.txt tail %r' % tail)
    return ''.join(lines[:-2])


def main():
    if shutil.which(T.MPIRUN) is None:
        print('FAIL: %s not found' % T.MPIRUN)
        return 1
    fails = []
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = T._make_case(tmp, 'reader', tension=False, caps=-1)
        path = os.path.join(case_dir, 'bGlobal.txt')
        base = _base_text(path)
        for label, tail, want in VARIANTS:
            with open(path, 'w') as f:
                f.write(base + tail)
            rc_f, out_f = T._run_fortran(case_dir)
            rc_p, out_p = T._run_python(case_dir)
            stale = want == ERR_INPUT_FILE_STALE
            bad = []
            if rc_f != want or (stale and STALE_TEXT not in out_f):
                bad.append(('Fortran', rc_f, out_f, 'exit %d' % want
                            + (' and %r' % STALE_TEXT if stale else '')))
            if stale:
                if rc_p == 0 or STALE_TEXT not in out_p:
                    bad.append(('python', rc_p, out_p, 'a failure saying %r' % STALE_TEXT))
            elif rc_p != want:
                bad.append(('python', rc_p, out_p, 'exit %d' % want))
            for tag, rc, out, exp in bad:
                fails.append('%s: %s exited %d, expected %s:\n%s'
                             % (label, tag, rc, exp, out[-1200:]))
            if not bad:
                print('  PASS  %s: Fortran exit %d, python %s'
                      % (label, rc_f, 'stale' if stale else 'exit %d' % rc_p))
    for f in fails:
        print('  FAIL ', f)
    print('%s: %d variants, %d failure(s)'
          % ('SUCCESS' if not fails else 'FAILURE', len(VARIANTS), len(fails)))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
