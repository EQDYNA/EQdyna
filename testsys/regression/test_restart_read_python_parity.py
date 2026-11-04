#! /usr/bin/env python3
"""
Parity guard (PROJECT_RULES.md rule 10; board row 16c) for the python-jax
restart read, eqdyna.netcdf_io.read_on_fault_restart, against the Fortran
reader netcdf_read_on_fault_eqdyna_restart.

Reuses the fixture builders of test_multifault_restart_eqquasi_layout.py (not
modified): the same EQquasi-layout files, the same expected() values, the same
two faults of different extents with one off-diagonal node each. For every
file the REAL Fortran reader (compiled driver) and the python reader are run
on the same bytes and must agree with expected() AND with each other:
  - ntotft = 2 EQquasi layout; ntotft = 1 rank-2 file (dz /= dx not refused);
  - slots 31-36 land in the vel master/slave fric slots and init_vel seeds v1
    from them at the right equations (the mode-2 velocity initialisation);
  - refusals, with the Fortran's exit code and message text: 1-fault file at
    ntotft = 2 (22), 3-fault file at ntotft = 2 (22), dz /= dx at ntotft = 2
    (19, dz > dx and dz < dx), missing file and missing variable (71).
Fails loudly (rule 2) if mpif90/netCDF-Fortran are unavailable.
"""
import os
import shutil
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'src', 'python'))
import testsys.regression.test_multifault_item7_fix as item7  # noqa: E402
import testsys.regression.test_multifault_restart_eqquasi_layout as lay  # noqa: E402
from eqdyna import assembleGlobalMass, netcdf_io  # noqa: E402
from eqdyna.checkInputConsistency import InputConsistencyError  # noqa: E402

NAME = 'test_restart_read_python_parity'
SLOTS = [8, 49, 7, 47, 20, 23, 31, 32, 33, 34, 35, 36]   # lay.VAR_NAMES order
# The driver's two nodes as the python reader sees them: 1-indexed node ids
# (meshCoor rows 1, 2), nsmp rows [slave, master, 0-indexed fault].
NSMP = np.array([[1, 2, 0], [2, 1, 1]], dtype=np.int64)


def py_read(nc, ntotft, dz):
    mesh = np.zeros((3, 3))
    for ift, f in enumerate(lay.FAULTS, start=1):
        mesh[ift, 0], mesh[ift, 2] = f['node']
    fric = np.full((3, 101), -999.0)
    boxes = [f['box'] for f in lay.FAULTS][:ntotft]
    return netcdf_io.read_on_fault_restart(
        nc, fric, ntotft,
        [b[0] for b in boxes], [b[1] for b in boxes], [b[2] for b in boxes], [b[3] for b in boxes],
        1.0, dz, mesh, NSMP[:ntotft])


def py_got(out, ntotft):
    return {ift: [float(out[ift, s]) for s in SLOTS] for ift in range(1, ntotft + 1)}


def msg(r):
    return r.stdout + r.stderr


def main():
    if shutil.which('mpif90') is None:
        print(f'FAIL {NAME}: mpif90 not on PATH (this guard requires a real Fortran compile)')
        return 1
    tmp = tempfile.mkdtemp(prefix='testRestartPyParity.')
    fails = []
    try:
        item7._DRIVER_SRC = lay.driver_src()
        try:
            binary = item7.build_driver(tmp)
        except RuntimeError as e:
            print(f'FAIL {NAME}\n -', e)
            return 1
        nc = os.path.join(tmp, 'fault.r.nc')
        exp = lay.expected()

        # values: ntotft = 2, EQquasi layout; Fortran == python == expected
        lay.write_eqquasi_layout(nc)
        r, fgot = lay.run(binary, tmp)
        if r.returncode != 0:
            fails.append(f'Fortran reader exited {r.returncode}: {msg(r)}')
        out2 = py_read(nc, 2, 1.0)
        pgot = py_got(out2, 2)
        for ift in (1, 2):
            for name, e, p in zip(lay.VAR_NAMES, exp[ift], pgot[ift]):
                if abs(p - e) > 1e-3:
                    fails.append(f'ntotft=2 python fault {ift} {name} = {p!r}, expected {e!r}')
            if len(fgot.get(ift, [])) != 12:
                fails.append(f'ntotft=2 fault {ift}: Fortran printed {fgot.get(ift)!r}')
                continue
            # Fortran prints ES25.16E3: the same float32 -> float64 value to 17 digits.
            for name, fv, p in zip(lay.VAR_NAMES, fgot[ift], pgot[ift]):
                if abs(fv - p) > 1e-13 * max(1.0, abs(p)):
                    fails.append(f'ntotft=2 fault {ift} {name}: Fortran {fv!r} != python {p!r}')
        touched = np.zeros(101, bool)
        touched[SLOTS] = True
        if not np.all(out2[1:, ~touched][:, 1:] == -999.0):
            fails.append('python reader wrote a slot outside the 12 restart slots')

        # init_vel seeds v1: slave from slots 34-36, master from 31-33
        eq_nums = {1: np.array([1, 2, 3]), 2: np.array([4, 5, 6])}
        v1 = assembleGlobalMass.init_vel(NSMP[:1, :2], eq_nums, out2[:2], 6)
        want = list(exp[1][9:12]) + list(exp[1][6:9])
        got = list(v1[1:7])
        if any(abs(a - b) > 1e-3 for a, b in zip(got, want)):
            fails.append(f'init_vel v1 = {got}, expected slave+master {want}')

        # ntotft = 1, rank-2 file, dz /= dx is NOT refused
        os.remove(nc)
        lay.write_rank2(nc)
        r, fgot = lay.run(binary, tmp, ('1', '2.0'))
        if r.returncode != 0:
            fails.append(f'Fortran rank-2 ntotft=1 dz=2 exited {r.returncode}: {msg(r)}')
        else:
            p1 = py_got(py_read(nc, 1, 2.0), 1)
            for name, fv, p in zip(lay.VAR_NAMES, fgot.get(1, []), p1[1]):
                if abs(fv - p) > 1e-13 * max(1.0, abs(p)):
                    fails.append(f'ntotft=1 dz=2 {name}: Fortran {fv!r} != python {p!r}')
            if len(fgot.get(1, [])) != 12:
                fails.append(f'ntotft=1 dz=2: Fortran printed {fgot!r}')

        def refusal(label, run_args, py_args, code, needles):
            r, _ = lay.run(binary, tmp, run_args)
            if r.returncode == 0:
                fails.append(f'{label}: Fortran returned 0')
            for n in needles:
                if n not in msg(r):
                    fails.append(f'{label}: Fortran output lacks {n!r}')
            try:
                py_read(*py_args)
                fails.append(f'{label}: python reader returned instead of raising')
            except InputConsistencyError as e:
                if e.code != code:
                    fails.append(f'{label}: python code {e.code}, expected {code}')
                for n in needles:
                    if n not in str(e):
                        fails.append(f'{label}: python message lacks {n!r}: {e}')

        os.remove(nc)
        lay.write_rank2(nc)
        refusal('1-fault file, ntotft=2', (), (nc, 2, 1.0), 22,
                ['fault.r.nc holds 1 fault(s) but ntotft = 2',
                 'write it with an EQquasi that carries the nid_fault dimension'])
        os.remove(nc)
        lay.write_eqquasi_layout(nc, extra_faults=1)
        refusal('3-fault file, ntotft=2', (), (nc, 2, 1.0), 22,
                ['fault.r.nc holds 3 fault(s) but ntotft = 2',
                 'it was written for a case with more faults than this one'])
        os.remove(nc)
        lay.write_eqquasi_layout(nc)
        for dz in (2.0, 0.5):
            refusal(f'dz={dz} ntotft=2', ('2', str(dz)), (nc, 2, dz), 19,
                    ['multi-fault restart (mode = 2) refused', 'dz = ', 'but dx = ',
                     'set dz equal to dx'])

        import netCDF4
        os.remove(nc)
        with netCDF4.Dataset(nc, 'w') as ds:
            ds.createDimension('d', 3)
            for name in lay.VAR_NAMES[:-1]:        # no vzs
                ds.createVariable(name, 'f4', ('d', 'd'))
        refusal('missing variable vzs', ('1', '1.0'), (nc, 1, 1.0), 71, ['NetCDF: '])
        os.remove(nc)
        refusal('missing file', ('1', '1.0'), (nc, 1, 1.0), 71, ['NetCDF: '])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if fails:
        print(f'FAIL {NAME}')
        for f in fails:
            print('  -', f)
        return 1
    print(f'SUCCESS {NAME} (python restart reader == Fortran reader on the EQquasi layout at '
          'ntotft=2 and the rank-2 file at ntotft=1, slots 31-36 seed init_vel; 1-fault, 3-fault, '
          'dz != dx, missing file, missing variable refused with the Fortran exit code and message)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
