#! /usr/bin/env python3
"""
Regression guard (issue #201): a rate-and-state run whose normal traction is
>= 0 at a fault node STOPS with ERR_NUM_NORMAL_STRESS_NONCOMPRESSIVE (65) and
prints the time, the node coordinates and the value. Before #201 both codes
clamped the value to 0 and carried on silently.

One tiny planar friclaw-4 case (dx=500 m, 3 steps, no nucleation), uniform
-120 MPa normal / 40 MPa shear, with ONE fault grid point set to +5 MPa
(tension) at x=500 m, z=-1000 m. Run through the REAL Fortran binary
(bin/eqdyna) and the REAL Python port (`python3 -m eqdyna ... --backend
numpy`):

  1. tension, explicit caps=-1 (off on a planar fault; the pre-board-commit-
     68fe99a default, still a valid explicit input, just no longer what a
     bare case gets)  -> both exit 65 and name the node (500, 0, -1000) and
     a value near +5 MPa;
  2. the same tension with par.C_normal_stress_caps = 1 -> both exit 0 (the
     cap lifts +5 MPa to min_norm = -10 MPa, so the switch reaches the
     solver and the stop sits AFTER the caps, as in the Fortran);
  3. a control with no tension node                     -> both exit 0;
  4. friclaw 5 (thermal pressurization), every node at -30 MPa normal and
     NO tension node: the pore pressure grows until the effective normal
     traction (tnrm + TP pore pressure, faulting.f90 solveRSF) turns
     positive -> both exit 65 at the same step, time, node and value, and
     the value is > 0 although the raw normal never is. Goes red if the
     stop moves above the TP offset (raw -30 MPa never stops) or if Python
     latches the pre-offset traction.

Slip-weakening (friclaw 1/2) is NOT stopped: solveSWTW keeps effNorm = 0
under tension (strength = cohesion), the SCEC treatment. That path is not
changed by #201 and every gated SW case pins it bitwise.

Fails (never skips) when mpirun or bin/eqdyna is missing.
"""
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpirun_capture  # noqa: E402  (rank-owned output past MPI_Abort)

MPIRUN = os.environ.get('EQDYNA_MPIRUN', 'mpirun')
CASE_NAME = 'test.tpv104'
ERR_NUM_NORMAL_STRESS_NONCOMPRESSIVE = 65
NSTEPS = 3
BAD_X, BAD_Z, BAD_TN = 500.0, -1000.0, 5.0e6

USER_PARAMS = '''#! /usr/bin/env python3
from defaultParameters import parameters
from math import *
from lib import *
import numpy as np

par = parameters()

par.xmin, par.xmax = -3.0e3, 3.0e3
par.ymin, par.ymax = -3.0e3, 3.1e3
par.zmin, par.zmax = -3.0e3, 0.0e3

par.fxmin, par.fxmax = -1.5e3, 1.5e3
par.fymin, par.fymax = 0.0, 0.0
par.fzmin, par.fzmax = -2.0e3, 0.0e3

par.xsource, par.ysource, par.zsource = 0.0, 0.0, -1000.0

par.dx = 500.
par.dy = par.dx
par.dz = par.dx
par.nmat = 1
par.vp, par.vs, par.rou = 6000., 3464., 2670.

par.dt = 0.5*par.dx/par.vp
par.term = %(nsteps)d.*par.dt
par.friclaw = 4
par.tpv = 104
par.C_nuclea = 0
par.C_normal_stress_caps = %(caps)d

par.nx, par.ny, par.nz = 1, 1, 1

par.nfx = round((par.fxmax - par.fxmin)/par.dx + 1)
par.nfz = round((par.fzmax - par.fzmin)/par.dz + 1)
par.fx  = np.linspace(par.fxmin,par.fxmax,par.nfx)
par.fz  = np.linspace(par.fzmin,par.fzmax,par.nfz)

par.on_fault_vars = np.zeros((par.nfz,par.nfx,100))
for ix, xcoor in enumerate(par.fx):
  for iz, zcoor in enumerate(par.fz):
    v = par.on_fault_vars[iz,ix]
    v[1], v[2], v[3] = par.fric_sw_fs, par.fric_sw_fd, par.fric_sw_D0
    v[7] = -120.0e6
    v[8] = 40.0e6
    v[9]  = par.fric_rsf_a
    v[10] = par.fric_rsf_b
    v[11] = par.fric_rsf_Dc
    v[12] = par.fric_rsf_v0
    v[13] = par.fric_rsf_r0
    v[14] = par.fric_rsf_fw
    v[15] = par.fric_rsf_vw
    v[46] = par.creep_slip_rate
    v[20] = state_steady_state(v[9], v[10], v[11], v[12], v[13], v[8], v[7],
                               v[46], par.friclaw)
    if %(tension)d and abs(xcoor-%(bx)r)<0.01 and abs(zcoor-%(bz)r)<0.01:
        v[7] = %(btn)r     # tension at one node; state kept from -120 MPa

par.st_coor_on_fault = [[0.0, -0.5]]
par.st_coor_off_fault = [[0,1,0]]
par.n_on_fault  = len(par.st_coor_on_fault)
par.n_off_fault = len(par.st_coor_off_fault)
'''


# Case 4: friclaw 5, uniform -30 MPa normal / 20 MPa shear, TP on, no
# tension node. TP_LAMBDA (on_fault_vars col 19) sets how fast the pore
# pressure rises; with it the effective normal turns positive at step 3.
TP_LAMBDA = 1.0e7
TP_NSTEPS = 5
TP_LOOP = """    v = par.on_fault_vars[iz,ix]
    v[1], v[2], v[3] = par.fric_sw_fs, par.fric_sw_fd, par.fric_sw_D0
    v[7] = -30.0e6
    v[8] = 20.0e6
    v[9]  = par.fric_rsf_a
    v[10] = par.fric_rsf_b
    v[11] = par.fric_rsf_Dc
    v[12] = par.fric_rsf_v0
    v[13] = par.fric_rsf_r0
    v[14] = par.fric_rsf_fw
    v[15] = par.fric_rsf_vw
    v[16] = par.fric_tp_a_hy
    v[17] = par.fric_tp_a_th
    v[18] = par.fric_tp_rouc
    v[19] = %r
    v[40] = par.fric_tp_h
    v[41] = par.fric_tp_Tini
    v[42] = par.fric_tp_pini
    v[46] = par.creep_slip_rate
    v[20] = state_steady_state(v[9], v[10], v[11], v[12], v[13], v[8], v[7],
                               1.0, par.friclaw)

""" % TP_LAMBDA


def _tp_params():
    """USER_PARAMS for case 4: friclaw 5 and the TP loop body."""
    up = USER_PARAMS.replace('par.friclaw = 4\npar.tpv = 104',
                             'par.friclaw = 5\npar.tpv = 105')
    i, j = up.index('    v = par.on_fault_vars[iz,ix]'), up.index('par.st_coor_on_fault')
    assert 'friclaw = 5' in up and j > i
    return up[:i] + TP_LOOP + up[j:]


def _env():
    env = dict(os.environ)
    env['EQDYNAROOT'] = ROOT
    env['PYTHONPATH'] = os.path.join(ROOT, 'src', 'python')
    env['PATH'] = os.pathsep.join([os.path.join(ROOT, 'bin'),
                                   os.path.join(ROOT, 'scripts'),
                                   env.get('PATH', '')])
    return env


def _make_case(tmp, name, tension, caps, nsteps=NSTEPS, extra='', params=None):
    """extra: user_defined_params lines appended after the caps switch.
    params: the user_defined_params template (default USER_PARAMS)."""
    case_dir = os.path.join(tmp, 'case_%s' % name)
    r = subprocess.run([sys.executable, os.path.join(ROOT, 'scripts', 'create.newcase'),
                        case_dir, CASE_NAME], env=_env(), capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError('create.newcase failed (%d):\n%s'
                           % (r.returncode, r.stdout[-2000:] + r.stderr[-2000:]))
    with open(os.path.join(case_dir, 'user_defined_params.py'), 'w') as f:
        f.write((params or USER_PARAMS).replace('par.C_normal_stress_caps = %(caps)d\n',
                                    'par.C_normal_stress_caps = %(caps)d\n' + extra)
                % dict(nsteps=nsteps, caps=caps, tension=int(tension),
                       bx=BAD_X, bz=BAD_Z, btn=BAD_TN))
    r = subprocess.run([sys.executable, 'case.setup'], cwd=case_dir,
                       env=_env(), capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError('case.setup failed for %s:\n%s' % (name, r.stdout + r.stderr))
    return case_dir


def _binary():
    for cand in (os.path.join(ROOT, 'bin', 'eqdyna'),
                 os.path.join(ROOT, 'src', 'fortran', 'eqdyna')):
        if os.path.exists(cand):
            return cand
    raise AssertionError('bin/eqdyna (or src/fortran/eqdyna) not found; '
                         'build it first with ./install-eqdyna.sh')


def _run_fortran(case_dir):
    return mpirun_capture.run_rank_files(MPIRUN, _binary(), case_dir, env=_env())


def _run_python(case_dir, nsteps=NSTEPS):
    r = subprocess.run([sys.executable, '-u', '-m', 'eqdyna', case_dir, str(nsteps),
                        '--backend', 'numpy'],
                       cwd=ROOT, env=_env(), capture_output=True, text=True, timeout=300)
    return r.returncode, r.stdout + r.stderr


NUM = r'([-+]?\d+\.\d+E[-+]\d+)'


def _check_report(tag, out):
    """The node coordinates and the value, from the text each side prints."""
    if tag == 'fortran':
        m = re.search(r'node x, y, z \(m\):\s*' + NUM + r'\s*' + NUM + r'\s*' + NUM, out)
        v = re.search(r'normal traction :\s*' + NUM + r' MPa', out)
        t = re.search(r'time \(s\)\s*:\s*' + NUM, out)
    else:
        m = re.search(r'node x, y, z \(m\) ' + NUM + ' ' + NUM + ' ' + NUM, out, re.I)
        v = re.search(r'normal traction ' + NUM + ' MPa', out, re.I)
        t = re.search(r'time ' + NUM + ' s', out, re.I)
    if not (m and v and t):
        raise AssertionError('%s: stop report lacks time/coordinates/value:\n%s'
                             % (tag, out[-2000:]))
    xyz = [float(g) for g in m.groups()]
    if max(abs(xyz[0] - BAD_X), abs(xyz[1]), abs(xyz[2] - BAD_Z)) > 1e-3:
        raise AssertionError('%s: stop names node %r, expected (%g, 0, %g)'
                             % (tag, xyz, BAD_X, BAD_Z))
    val = float(v.group(1))
    if not (0.0 <= val and abs(val - BAD_TN / 1e6) < 1.0):
        raise AssertionError('%s: stop value %g MPa, expected ~ +%g'
                             % (tag, val, BAD_TN / 1e6))
    if float(t.group(1)) <= 0.0:
        raise AssertionError('%s: stop time %s is not > 0' % (tag, t.group(1)))
    return xyz, val


def check_stop_both_codes():
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = _make_case(tmp, 'tension', tension=True, caps=-1)
        rc_f, out_f = _run_fortran(case_dir)
        if rc_f != ERR_NUM_NORMAL_STRESS_NONCOMPRESSIVE:
            raise AssertionError('Fortran exited %d on a +5 MPa RSF node, expected %d:\n%s'
                                 % (rc_f, ERR_NUM_NORMAL_STRESS_NONCOMPRESSIVE, out_f[-2000:]))
        if 'Non-compressive normal traction' not in out_f:
            raise AssertionError('Fortran exit 65 without the diagnostic:\n%s' % out_f[-2000:])
        xf, vf = _check_report('fortran', out_f)
        rc_p, out_p = _run_python(case_dir)
        if rc_p != ERR_NUM_NORMAL_STRESS_NONCOMPRESSIVE:
            raise AssertionError('python exited %d on a +5 MPa RSF node, expected %d:\n%s'
                                 % (rc_p, ERR_NUM_NORMAL_STRESS_NONCOMPRESSIVE, out_p[-2000:]))
        if 'exit code :  65' not in out_p:
            raise AssertionError('python exit 65 without the FATAL block:\n%s' % out_p[-2000:])
        xp_, vp = _check_report('python', out_p)
        if abs(vf - vp) > 1e-6:
            raise AssertionError('Fortran value %.8e MPa != python %.8e MPa' % (vf, vp))
    print('  PASS  +5 MPa RSF node: Fortran and python both exit 65 at node '
          '(%g, %g, %g), %.4f MPa' % (xf[0], xf[1], xf[2], vf))


def _parse_report(tag, out):
    """(friclaw, time, step, (x, y, z), value MPa) from either code's report."""
    if tag == 'fortran':
        fl = re.search(r'friclaw (\d+)', out)
        ts = re.search(r'time \(s\)\s*:\s*' + NUM + r'\s+step\s+(\d+)', out)
        m = re.search(r'node x, y, z \(m\):\s*' + NUM + r'\s*' + NUM + r'\s*' + NUM, out)
        v = re.search(r'normal traction :\s*' + NUM + r' MPa', out)
    else:
        fl = re.search(r'friclaw (\d+)', out)
        ts = re.search(r'time ' + NUM + r' s, step (\d+)', out, re.I)
        m = re.search(r'node x, y, z \(m\) ' + NUM + ' ' + NUM + ' ' + NUM, out, re.I)
        v = re.search(r'normal traction ' + NUM + ' MPa', out, re.I)
    if not (fl and ts and m and v):
        raise AssertionError('%s: stop report lacks friclaw/time/step/coordinates/'
                             'value:\n%s' % (tag, out[-2000:]))
    return (int(fl.group(1)), float(ts.group(1)), int(ts.group(2)),
            tuple(float(g) for g in m.groups()), float(v.group(1)))


def check_tp_stop_both_codes():
    """Case 4: the stop sees the TP-offset traction, in both codes."""
    with tempfile.TemporaryDirectory() as tmp:
        case_dir = _make_case(tmp, 'tp', tension=False, caps=-1,
                              nsteps=TP_NSTEPS, params=_tp_params())
        reps = {}
        for tag, run in (('fortran', lambda: _run_fortran(case_dir)),
                         ('python', lambda: _run_python(case_dir, nsteps=TP_NSTEPS))):
            rc, out = run()
            if rc != ERR_NUM_NORMAL_STRESS_NONCOMPRESSIVE:
                raise AssertionError(
                    '%s exited %d on friclaw 5 with the pore pressure past a -30 MPa '
                    'normal, expected %d (does the stop see the TP offset?):\n%s'
                    % (tag, rc, ERR_NUM_NORMAL_STRESS_NONCOMPRESSIVE, out[-2000:]))
            reps[tag] = _parse_report(tag, out)
    (lf, tf, sf, xf, vf), (lp, tp, sp, xp_, vp) = reps['fortran'], reps['python']
    if lf != 5 or lp != 5:
        raise AssertionError('friclaw printed %d (Fortran) / %d (python), expected 5'
                             % (lf, lp))
    if not vf > 0.0:
        raise AssertionError('Fortran stop value %g MPa is not > 0: the raw normal '
                             'is -30 MPa, only the TP offset makes it >= 0' % vf)
    if sf != sp or abs(tf - tp) > 1e-12 * max(1.0, abs(tf)):
        raise AssertionError('stop at step %d, t=%.8e (Fortran) vs step %d, t=%.8e '
                             '(python)' % (sf, tf, sp, tp))
    if max(abs(a - b) for a, b in zip(xf, xp_)) > 1e-3:
        raise AssertionError('stop node %r (Fortran) vs %r (python)' % (xf, xp_))
    if abs(vf - vp) > 1e-6 * max(1.0, abs(vf)):
        raise AssertionError('Fortran value %.8e MPa != python %.8e MPa' % (vf, vp))
    print('  PASS  friclaw 5 TP: both exit 65 at step %d, t=%.4e s, node '
          '(%g, %g, %g), %+.4f MPa (raw normal -30 MPa)'
          % (sf, tf, xf[0], xf[1], xf[2], vf))


def check_caps_and_control_run():
    for name, tension, caps in (('caps_on', True, 1), ('control', False, -1)):
        with tempfile.TemporaryDirectory() as tmp:
            case_dir = _make_case(tmp, name, tension=tension, caps=caps)
            rc_f, out_f = _run_fortran(case_dir)
            if rc_f != 0:
                raise AssertionError('%s: Fortran exited %d, expected 0:\n%s'
                                     % (name, rc_f, out_f[-2000:]))
            rc_p, out_p = _run_python(case_dir)
            if rc_p != 0:
                raise AssertionError('%s: python exited %d, expected 0:\n%s'
                                     % (name, rc_p, out_p[-2000:]))
        print('  PASS  %s: Fortran and python both exit 0' % name)


def main():
    import shutil
    if shutil.which(MPIRUN) is None:
        print('FAIL: %s not found' % MPIRUN)
        return 1
    fails = 0
    for check in (check_stop_both_codes, check_caps_and_control_run,
                  check_tp_stop_both_codes):
        try:
            check()
        except Exception as e:  # noqa: BLE001 -- report every check
            fails += 1
            print('  FAIL  %s: %s' % (check.__name__, e))
    print('%s: %d failure(s)' % ('SUCCESS' if fails == 0 else 'FAILURE', fails))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
