#! /usr/bin/env python3
"""
Regression guard (PROJECT_RULES.md rule 10) for the multi-fault restart
handoff between EQquasi and EQdyna in a fully dynamic cycle (mode 2).

THE BUG: EQquasi's netcdf_write_on_fault writes fault.r.nc as ONE set of
untagged variables with an explicit nid_fault dimension -- Fortran dims
(nid_strike, nid_dip, nid_fault), each fault's own (strike, dip) grid at the
(1,1) corner of its slice, zero-padded to the largest fault. EQdyna's
netcdf_read_on_fault_eqdyna_restart instead looked up faultTag()-prefixed
2-D variables ('ft2_shear_strike', ...) that no code writes, so any
ntotft>1 restart died in nf90_inq_varid. In the other direction
plotRuptureDynamics wrote fault 2 to a separate fault.dyna.r_ft2.nc, which
EQquasi's netcdf_read_on_fault_restart never opens: it reads fault.dyna.r.nc
with the same nid_fault layout it writes.

THIS TEST compiles a tiny driver against the REAL netcdf_io.f90 (same build
as test_multifault_item7_fix.py), places two faults of DIFFERENT extents
(fault 1: 3 strike x 2 dip nodes; fault 2: 2 x 3, so the file is padded to
3 x 3) with one node each at an OFF-DIAGONAL cell, and checks the fric slots
the real reader fills. Every cell holds a distinct value encoding
(fault, variable, strike, dip), so a transpose, a wrong slice, or one fault
aliased onto the other all read a different number. Three files go in:
  1. fault.r.nc built in EQquasi's layout (what EQquasi hands EQdyna);
  2. fault.dyna.r.nc from the REAL scripts/plotRuptureDynamics
     generateNcRestart, read back through the same reader -- the writer must
     produce the nid_fault layout, not a second per-fault file;
  3. a rank-2 (pre-nid_fault) file with ntotft=2, which must abort non-zero
     rather than give fault 2 fault 1's data.

Cheap (rule 9): four small .f90 compiles and three tiny netCDF reads, a few
seconds. Fails loudly (rule 2) if mpif90/netCDF-Fortran are unavailable.
"""
import importlib.machinery
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
import types

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
import testsys.regression.test_multifault_item7_fix as item7  # noqa: E402  (same build)

NAME = 'test_multifault_restart_eqquasi_layout'

# The 12 restart variables, in netcdf_read_on_fault_eqdyna_restart's order.
VAR_NAMES = ['shear_strike', 'shear_dip', 'effective_normal', 'slip_rate',
             'state_variable', 'state_normal', 'vxm', 'vym', 'vzm',
             'vxs', 'vys', 'vzs']

# Per fault: (fxmin, fxmax, fzmin, fzmax), dx = dz = 1, and the one node's
# (x, z), chosen off the diagonal of the fault's own grid.
FAULTS = [
    dict(box=(0.0, 2.0, -1.0, 0.0), node=(2.0, -1.0)),   # 3x2, node (ii,jj)=(3,1)
    dict(box=(5.0, 6.0, -2.0, 0.0), node=(6.0, -1.0)),   # 2x3, node (ii,jj)=(2,2)
]

_DRIVER_SRC = r'''
program restart_layout_driver
    use globalvar
    implicit none
    integer :: ift, k

    ntotft = 2
    dx = 1.0d0
    dz = 1.0d0

    allocate(fxmin(2), fxmax(2), fzmin(2), fzmax(2))
%(boxes)s
    allocate(nftnd(2))
    nftnd = 1
    allocate(nsmp(2,1,2))
    nsmp(1,1,1) = 1
    nsmp(1,1,2) = 2
    allocate(meshCoor(3,2))
%(nodes)s
    allocate(fric(100,1,2))
    fric = -999.0d0

    call netcdf_read_on_fault_eqdyna_restart

    do ift = 1, 2
        write(*,'(A,I0,1X,12ES25.16E3)') 'FAULT', ift, &
            fric(FRIC_SLOT_INIT_STRIKE_SHEAR, 1, ift), fric(FRIC_SLOT_INIT_DIP_SHEAR, 1, ift), &
            fric(FRIC_SLOT_INIT_NORM, 1, ift), fric(FRIC_SLOT_PEAK_SLIPRATE, 1, ift), &
            fric(FRIC_SLOT_STATE, 1, ift), fric(FRIC_SLOT_THETA_PC, 1, ift), &
            fric(FRIC_SLOT_VEL_MASTER_X, 1, ift), fric(FRIC_SLOT_VEL_MASTER_Y, 1, ift), &
            fric(FRIC_SLOT_VEL_MASTER_Z, 1, ift), fric(FRIC_SLOT_VEL_SLAVE_X, 1, ift), &
            fric(FRIC_SLOT_VEL_SLAVE_Y, 1, ift), fric(FRIC_SLOT_VEL_SLAVE_Z, 1, ift)
    enddo
end program restart_layout_driver
'''


def driver_src():
    boxes, nodes = [], []
    for ift, f in enumerate(FAULTS, start=1):
        xlo, xhi, zlo, zhi = f['box']
        boxes.append(f'    fxmin({ift}) = {xlo}d0; fxmax({ift}) = {xhi}d0; '
                     f'fzmin({ift}) = {zlo}d0; fzmax({ift}) = {zhi}d0')
        nodes.append(f'    meshCoor(1,{ift}) = {f["node"][0]}d0; '
                     f'meshCoor(3,{ift}) = {f["node"][1]}d0')
    return _DRIVER_SRC % dict(boxes='\n'.join(boxes), nodes='\n'.join(nodes))


def grid(f):
    xlo, xhi, zlo, zhi = f['box']
    return int(round(xhi - xlo)) + 1, int(round(zhi - zlo)) + 1   # nfx, nfz


def value(ift, k, i, j):
    """Distinct per (fault 1.., variable 1.., strike 1.., dip 1..); exact in float32."""
    return 1000.0*ift + 10.0*k + 0.1*i + 0.01*j


def fault_arrays():
    """Each fault's (nfz, nfx, 12) array, the (dip, strike, var) shape
    loadFrtData's fVarArr and EQquasi's per-fault slab both use."""
    out = []
    for ift, f in enumerate(FAULTS, start=1):
        nfx, nfz = grid(f)
        a = np.zeros((nfz, nfx, len(VAR_NAMES)))
        for k in range(len(VAR_NAMES)):
            for i in range(nfx):
                for j in range(nfz):
                    a[j, i, k] = value(ift, k + 1, i + 1, j + 1)
        out.append(a)
    return out


def write_eqquasi_layout(path):
    """fault.r.nc as EQquasi's netcdf_write_on_fault writes it: NF90_REAL,
    C-order dims (nid_fault, nid_dip, nid_strike), zero padding."""
    import netCDF4
    arrs = fault_arrays()
    nfxMax = max(a.shape[1] for a in arrs)
    nfzMax = max(a.shape[0] for a in arrs)
    with netCDF4.Dataset(path, 'w') as ds:
        ds.createDimension('nid_dip', nfzMax)
        ds.createDimension('nid_strike', nfxMax)
        ds.createDimension('nid_fault', len(arrs))
        for k, name in enumerate(VAR_NAMES):
            v = ds.createVariable(name, 'f4', ('nid_fault', 'nid_dip', 'nid_strike'))
            data = np.zeros((len(arrs), nfzMax, nfxMax))
            for ift, a in enumerate(arrs):
                data[ift, :a.shape[0], :a.shape[1]] = a[:, :, k]
            v[:] = data


def write_rank2(path):
    """A pre-nid_fault file: one fault's 2-D variables only."""
    import netCDF4
    a = fault_arrays()[0]
    with netCDF4.Dataset(path, 'w') as ds:
        ds.createDimension('dip', a.shape[0])
        ds.createDimension('strike', a.shape[1])
        for k, name in enumerate(VAR_NAMES):
            ds.createVariable(name, 'f8', ('dip', 'strike'))[:] = a[:, :, k]


def write_with_plotRuptureDynamics(tmp):
    """fault.dyna.r.nc from the real generateNcRestart. The script imports
    user_defined_params.par and lib at module level; par is not used by
    generateNcRestart, so a stub module is enough."""
    sys.modules.setdefault('user_defined_params',
                           types.SimpleNamespace(par=types.SimpleNamespace()))
    sys.path.insert(0, os.path.join(ROOT, 'scripts'))
    import matplotlib
    matplotlib.use('Agg')
    path = os.path.join(ROOT, 'scripts', 'plotRuptureDynamics')
    loader = importlib.machinery.SourceFileLoader('plotRuptureDynamics', path)
    spec = importlib.util.spec_from_loader('plotRuptureDynamics', loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    faults = []
    for f, a in zip(FAULTS, fault_arrays()):
        nfx, nfz = grid(f)
        padded = np.zeros((nfz, nfx, 100))
        padded[:, :, :len(VAR_NAMES)] = a
        faults.append((padded, None, None, nfx, nfz, ''))
    cwd = os.getcwd()
    os.chdir(tmp)
    try:
        mod.generateNcRestart(faults)
    finally:
        os.chdir(cwd)
    names = sorted(n for n in os.listdir(tmp) if n.startswith('fault.dyna.r'))
    if names != ['fault.dyna.r.nc']:
        raise RuntimeError(f'generateNcRestart wrote {names}, expected one fault.dyna.r.nc')
    os.replace(os.path.join(tmp, 'fault.dyna.r.nc'), os.path.join(tmp, 'fault.r.nc'))


def expected():
    out = {}
    for ift, f in enumerate(FAULTS, start=1):
        xlo, _, zlo, _ = f['box']
        i = int(round(f['node'][0] - xlo)) + 1
        j = int(round(f['node'][1] - zlo)) + 1
        out[ift] = [value(ift, k, i, j) for k in range(1, len(VAR_NAMES) + 1)]
    return out


def run(binary, tmp):
    r = subprocess.run([binary], cwd=tmp, capture_output=True, text=True, timeout=30)
    got = {}
    for line in r.stdout.splitlines():
        parts = line.split()
        if parts and parts[0].startswith('FAULT'):
            got[int(parts[0][5:])] = [float(p) for p in parts[1:]]
    return r, got


def check(label, got, fails):
    for ift, exp in expected().items():
        g = got.get(ift)
        if g is None or len(g) != len(exp):
            fails.append(f'{label}: fault {ift} missing from driver output ({got!r})')
            continue
        for name, e, v in zip(VAR_NAMES, exp, g):
            if abs(v - e) > 1e-3:
                fails.append(f'{label}: fault {ift} {name} = {v!r}, expected {e!r}')


def main():
    if shutil.which('mpif90') is None:
        print(f'FAIL {NAME}: mpif90 not on PATH '
              '(this guard requires a real Fortran compile, it does not skip silently)')
        return 1

    tmp = tempfile.mkdtemp(prefix='testRestartLayout.')
    fails = []
    try:
        item7._DRIVER_SRC = driver_src()   # item7.build_driver compiles this source
        try:
            binary = item7.build_driver(tmp)
        except RuntimeError as e:
            print(f'FAIL {NAME}')
            print(' -', e)
            return 1
        nc = os.path.join(tmp, 'fault.r.nc')

        write_eqquasi_layout(nc)
        r, got = run(binary, tmp)
        if r.returncode != 0:
            fails.append(f'EQquasi layout: driver exited {r.returncode}:\n{r.stdout}\n{r.stderr}')
        else:
            check('EQquasi layout', got, fails)

        os.remove(nc)
        write_with_plotRuptureDynamics(tmp)
        r, got = run(binary, tmp)
        if r.returncode != 0:
            fails.append(f'generateNcRestart output: driver exited {r.returncode}:\n{r.stdout}\n{r.stderr}')
        else:
            check('generateNcRestart output', got, fails)

        os.remove(nc)
        write_rank2(nc)
        r, got = run(binary, tmp)
        if r.returncode == 0:
            fails.append('rank-2 file with ntotft=2: reader returned 0 '
                         f'instead of aborting (fric: {got!r})')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if fails:
        print(f'FAIL {NAME}')
        for f in fails:
            print('  -', f)
        return 1
    print(f'SUCCESS {NAME} (EQquasi nid_fault layout read per fault; '
          'generateNcRestart writes the same layout; rank-2 multi-fault file refused)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
