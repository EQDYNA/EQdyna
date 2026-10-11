#! /usr/bin/env python3
"""
Two-fault restart round trip, EQdyna -> EQquasi -> EQdyna (PR #172).

  1. build THIS checkout's src/fortran into <workdir>/build (never into the tree);
  2. run test.tpv22 (two faults) for --steps time steps at --ranks MPI ranks,
     then the real scripts/plotRuptureDynamics -> fault.dyna.r.nc;
  3. copy fault.dyna.r.nc to fault.r.nc, exactly as EQquasi's
     scripts/case.setup does for cycle i > 1;
  4. a driver compiled against EQquasi's REAL src/globalvar.f90 and
     src/netcdf_io.f90 (compiled from --eqquasi-src, objects in <workdir>,
     nothing written to that tree) calls netcdf_read_on_fault_restart, then
     netcdf_write_on_fault -> a new fault.r.nc (NF90_REAL);
  5. a driver compiled against THIS checkout's netcdf_io.f90 (the same build as
     testsys/regression/test_multifault_restart_eqquasi_layout.py) calls
     netcdf_read_on_fault_eqdyna_restart on that file;
  6. every fault node in the run's frt.txt* is checked: the 12 restart values
     EQdyna reads back must equal float32(frt value) EXACTLY, per fault.

Stand-ins (disclosed, not hidden):
  * EQquasi takes zero quasi-static steps. Its reader fills the *0 / VINIT
    slots and its writer reads the current-state slots, so the driver copies
    FR_TSTK0->FR_TSTK, FR_TDIP0->FR_TDIP, FR_TNRM0->FR_TNRM, FR_VINIT->FR_V_TRIAL
    between the two calls (the slots a zero-length EQquasi run would leave
    equal). Every other restart slot is the same slot on both sides.
  * EQquasi's infile1 (a, b, Dc, v0, r0) is a synthetic file in EQquasi's own
    3-D layout; those five variables are not part of the handoff under test.
  * The node set is the fault nodes of the actual run (frt.txt* rows), each
    fault's (x, z) mapped by the readers' own index formulas. tpv22 has
    dx == dz, so EQquasi's dx-only dip index (see PR #172) does not bite here.

Usage:
  roundtrip.py --eqquasi-src <eqquasi>/src --workdir <scratch dir> [--steps 20 --ranks 4]
Prints the report to stdout; exit 0 only if every value matches.
"""
import argparse
import glob
import hashlib
import os
import shutil
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)
import testsys.regression.test_multifault_item7_fix as item7  # noqa: E402  (EQdyna netcdf_io build)

CASE = 'test.tpv22'
VARS = ['shear_strike', 'shear_dip', 'effective_normal', 'slip_rate',
        'state_variable', 'state_normal', 'vxm', 'vym', 'vzm', 'vxs', 'vys', 'vzs']
# frt.txt column of each restart variable (scripts/lib.py loadFrtData, fVarArr).
FRT_COL = [12, 13, 11, 10, 20, 21, 14, 15, 16, 17, 18, 19]

EQQ_DRIVER = r'''
program eqquasi_roundtrip
    use globalvar
    implicit none
    integer :: ift, n, k, nmax, ntot
    character (len = 50) :: f1, f2, fo
    open(10, file='nodes.txt', status='old')
    read(10,*) ntotft, dx
    allocate(fltxyz(2,4,ntotft), nftnd(ntotft))
    fltxyz = 0.0d0
    do ift = 1, ntotft
        read(10,*) fltxyz(1,1,ift), fltxyz(2,1,ift), fltxyz(1,3,ift), fltxyz(2,3,ift), nftnd(ift)
    enddo
    nmax = maxval(nftnd)
    ntot = sum(nftnd)
    allocate(nsmp(2,nmax,ntotft), x(3,ntot), fric(100,nmax,ntotft))
    nsmp = 0
    x = 0.0d0
    fric = 0.0d0
    k = 0
    do ift = 1, ntotft
        do n = 1, nftnd(ift)
            k = k + 1
            read(10,*) x(1,k), x(3,k)
            nsmp(1,n,ift) = k
            nsmp(2,n,ift) = k
        enddo
    enddo
    close(10)
    f1 = 'on_fault_vars_input.nc'
    f2 = 'fault.r.nc'
    call netcdf_read_on_fault_restart(f1, f2)
    do ift = 1, ntotft
        do n = 1, nftnd(ift)
            fric(FR_TSTK, n, ift)    = fric(FR_TSTK0, n, ift)
            fric(FR_TDIP, n, ift)    = fric(FR_TDIP0, n, ift)
            fric(FR_TNRM, n, ift)    = fric(FR_TNRM0, n, ift)
            fric(FR_V_TRIAL, n, ift) = fric(FR_VINIT, n, ift)
        enddo
    enddo
    fo = 'out/fault.r.nc'
    call netcdf_write_on_fault(fo)
    print *, 'EQQUASI_OK', ntotft, ntot
end program eqquasi_roundtrip
'''

EQD_DRIVER = r'''
program eqdyna_roundtrip
    use globalvar
    implicit none
    integer :: ift, n, k, nmax, ntot
    open(10, file='nodes.txt', status='old')
    read(10,*) ntotft, dx, dz
    allocate(fxmin(ntotft), fxmax(ntotft), fzmin(ntotft), fzmax(ntotft), nftnd(ntotft))
    do ift = 1, ntotft
        read(10,*) fxmin(ift), fxmax(ift), fzmin(ift), fzmax(ift), nftnd(ift)
    enddo
    nmax = maxval(nftnd)
    ntot = sum(nftnd)
    allocate(nsmp(2,nmax,ntotft), meshCoor(3,ntot), fric(100,nmax,ntotft))
    nsmp = 0
    meshCoor = 0.0d0
    fric = -999.0d0
    k = 0
    do ift = 1, ntotft
        do n = 1, nftnd(ift)
            k = k + 1
            read(10,*) meshCoor(1,k), meshCoor(3,k)
            nsmp(1,n,ift) = k
            nsmp(2,n,ift) = k
        enddo
    enddo
    close(10)
    call netcdf_read_on_fault_eqdyna_restart
    open(20, file='fric_out.txt', status='replace')
    do ift = 1, ntotft
        do n = 1, nftnd(ift)
            write(20,'(I0,1X,I0,12ES25.16E3)') ift, n, &
                fric(FRIC_SLOT_INIT_STRIKE_SHEAR, n, ift), fric(FRIC_SLOT_INIT_DIP_SHEAR, n, ift), &
                fric(FRIC_SLOT_INIT_NORM, n, ift), fric(FRIC_SLOT_PEAK_SLIPRATE, n, ift), &
                fric(FRIC_SLOT_STATE, n, ift), fric(FRIC_SLOT_THETA_PC, n, ift), &
                fric(FRIC_SLOT_VEL_MASTER_X, n, ift), fric(FRIC_SLOT_VEL_MASTER_Y, n, ift), &
                fric(FRIC_SLOT_VEL_MASTER_Z, n, ift), fric(FRIC_SLOT_VEL_SLAVE_X, n, ift), &
                fric(FRIC_SLOT_VEL_SLAVE_Y, n, ift), fric(FRIC_SLOT_VEL_SLAVE_Z, n, ift)
        enddo
    enddo
    close(20)
    print *, 'EQDYNA_OK', ntotft, ntot
end program eqdyna_roundtrip
'''


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def run(cmd, cwd, env=None, log=None):
    t0 = time.time()
    r = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True)
    if log:
        with open(os.path.join(cwd, log), 'w') as f:
            f.write(r.stdout + r.stderr)
    if r.returncode != 0:
        sys.stdout.write(r.stdout[-4000:] + r.stderr[-4000:])
        raise SystemExit('FAIL: `%s` exited %d' % (' '.join(os.path.basename(c) for c in cmd), r.returncode))
    return r, time.time() - t0


def git_head(path):
    return subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=path, capture_output=True,
                          text=True).stdout.strip()


def git_dirty(path):
    return bool(subprocess.run(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=path,
                               capture_output=True, text=True).stdout.strip())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--eqquasi-src', required=True)
    ap.add_argument('--workdir', required=True)
    ap.add_argument('--steps', type=int, default=20)
    ap.add_argument('--ranks', type=int, default=4)
    a = ap.parse_args()
    eqq = os.path.abspath(a.eqquasi_src)
    wd = os.path.abspath(a.workdir)
    if os.path.exists(wd):
        shutil.rmtree(wd)
    os.makedirs(wd)
    import netCDF4

    print('== inputs')
    print('EQdyna HEAD     ', git_head(ROOT), '(dirty tracked files)' if git_dirty(ROOT) else '(clean)')
    print('EQquasi HEAD    ', git_head(eqq), '(dirty tracked files)' if git_dirty(eqq) else '(clean)')
    for label, p in [('EQdyna src/fortran/netcdf_io.f90', os.path.join(ROOT, 'src/fortran/netcdf_io.f90')),
                     ('EQdyna scripts/plotRuptureDynamics', os.path.join(ROOT, 'scripts/plotRuptureDynamics')),
                     ('EQdyna scripts/lib.py', os.path.join(ROOT, 'scripts/lib.py')),
                     ('EQdyna case_input/%s/user_defined_params.py' % CASE,
                      os.path.join(ROOT, 'case_input', CASE, 'user_defined_params.py')),
                     ('EQdyna case_input/%s/tpv22_23_common.py' % CASE,
                      os.path.join(ROOT, 'case_input', CASE, 'tpv22_23_common.py')),
                     ('EQquasi src/netcdf_io.f90', os.path.join(eqq, 'netcdf_io.f90')),
                     ('EQquasi src/globalvar.f90', os.path.join(eqq, 'globalvar.f90'))]:
        print('sha256 %s  %s' % (sha(p), label))

    # 1. build eqdyna from this checkout, outside the tree.
    # The makefile stamps the binary with ../../scripts/src_hash.py, so the
    # copy keeps the tree's src/fortran + scripts/src_hash.py layout.
    build = os.path.join(wd, 'tree', 'src', 'fortran')
    shutil.copytree(os.path.join(ROOT, 'src', 'fortran'), build,
                    ignore=shutil.ignore_patterns('*.o', '*.mod', 'eqdyna', 'srcStamp.inc'))
    os.makedirs(os.path.join(wd, 'tree', 'scripts'))
    shutil.copy(os.path.join(ROOT, 'scripts', 'src_hash.py'), os.path.join(wd, 'tree', 'scripts'))
    env = dict(os.environ, MACHINE='ubuntu')
    _, tb = run(['make'], build, env, log='make.log')
    eqdyna = os.path.join(build, 'eqdyna')
    print('== 1. build: %.0f s, sha256 %s  eqdyna' % (tb, sha(eqdyna)))

    # 2. tpv22 for --steps steps.
    env = dict(os.environ, EQDYNAROOT=ROOT,
               PATH=os.path.join(ROOT, 'scripts') + os.pathsep + os.environ['PATH'])
    run(['create.newcase', 'case', CASE], wd, env)
    case = os.path.join(wd, 'case')
    with open(os.path.join(case, 'user_defined_params.py'), 'a') as f:
        f.write('\n# roundtrip.py: short run, %d steps\npar.term = %d * par.dt\n' % (a.steps, a.steps))
    run(['./case.setup'], case, env, log='case.setup.log')
    _, tr = run(['mpirun', '-np', str(a.ranks), eqdyna], case, env, log='eqdyna.log')
    _, tp = run([sys.executable, 'plotRuptureDynamics', '--skip-final-surf-disp'], case, env,
                log='plotRuptureDynamics.log')
    sys.path.insert(0, case)
    from user_defined_params import par
    geom = par.faultgeom
    ntotft = len(geom)
    nstep = None
    for line in open(os.path.join(case, 'eqdyna.log')):
        if 'nstep' in line.lower() or 'time step' in line.lower():
            nstep = line.strip()
    print('== 2. %s: ntotft=%d dx=%g dz=%g term=%.6g s (%d x dt=%.6g), %d ranks, solve %.0f s, '
          'plotRuptureDynamics %.0f s' % (CASE, ntotft, par.dx, par.dz, par.term, a.steps, par.dt,
                                          a.ranks, tr, tp))
    if nstep:
        print('   eqdyna log:', nstep)

    # frt rows -> per-fault node set and the original values (last row wins, as loadFrtData).
    frt = sorted(glob.glob(os.path.join(case, 'frt.txt*')))
    nodes = [dict() for _ in range(ntotft)]
    for p in frt:
        d = np.loadtxt(p, ndmin=2)
        for row in d:
            hit = [k for k, g in enumerate(geom) if g[2] - 1.0 <= row[1] <= g[3] + 1.0]
            if len(hit) != 1:
                raise SystemExit('FAIL: frt row at y=%r matches %d faults' % (row[1], len(hit)))
            nodes[hit[0]][(round(row[0], 3), round(row[2], 3))] = row[FRT_COL]
    print('   frt files: %d; fault nodes per fault: %s' % (len(frt), [len(n) for n in nodes]))
    dyna = os.path.join(case, 'fault.dyna.r.nc')
    with netCDF4.Dataset(dyna) as ds:
        print('   fault.dyna.r.nc dims %s, sha256 %s' % (
            {k: len(v) for k, v in ds.dimensions.items()}, sha(dyna)))

    def write_nodes(path, with_dz):
        with open(path, 'w') as f:
            # %.17g: exact float64 round trip (numpy 2's %r gives 'np.float64(..)')
            f.write('%d %.17g%s\n' % (ntotft, par.dx, (' %.17g' % par.dz) if with_dz else ''))
            for g, n in zip(geom, nodes):
                f.write('%.17g %.17g %.17g %.17g %d\n' % (g[0], g[1], g[4], g[5], len(n)))
            for n in nodes:
                for (x, z) in n:
                    f.write('%.17g %.17g\n' % (x, z))

    # 3. EQquasi stage: fault.dyna.r.nc -> fault.r.nc (case.setup:370-371), synthetic infile1.
    qd = os.path.join(wd, 'eqquasi_stage')
    os.makedirs(os.path.join(qd, 'out'))
    shutil.copy(dyna, os.path.join(qd, 'fault.r.nc'))
    nfx = [int(round((g[1] - g[0]) / par.dx)) + 1 for g in geom]
    nfz = [int(round((g[5] - g[4]) / par.dx)) + 1 for g in geom]   # EQquasi's own dx-only formula
    with netCDF4.Dataset(os.path.join(qd, 'on_fault_vars_input.nc'), 'w') as ds:
        ds.createDimension('nid_fault', ntotft)
        ds.createDimension('nid_dip', max(nfz))
        ds.createDimension('nid_strike', max(nfx))
        for k, name in enumerate(['a', 'b', 'Dc', 'v0', 'r0']):
            ds.createVariable(name, 'f4', ('nid_fault', 'nid_dip', 'nid_strike'))[:] = 0.01 * (k + 1)
    write_nodes(os.path.join(qd, 'nodes.txt'), with_dz=False)

    # 4. EQquasi driver against EQquasi's real globalvar.f90 + netcdf_io.f90.
    ob = os.path.join(wd, 'eqquasi_build')
    os.makedirs(ob)
    inc = subprocess.run(['nf-config', '--fflags'], capture_output=True, text=True).stdout.split()
    lib = subprocess.run(['nf-config', '--flibs'], capture_output=True, text=True).stdout.split()
    fc = ['mpif90', '-O0', '-ffree-line-length-none', '-I', ob, '-J', ob]
    run(fc + ['-c', os.path.join(eqq, 'globalvar.f90'), '-o', 'globalvar.o'], ob)
    run(fc + inc + ['-c', os.path.join(eqq, 'netcdf_io.f90'), '-o', 'netcdf_io.o'], ob)
    with open(os.path.join(ob, 'driver.f90'), 'w') as f:
        f.write(EQQ_DRIVER)
    run(fc + ['-c', 'driver.f90', '-o', 'driver.o'], ob)
    run(['mpif90', 'globalvar.o', 'netcdf_io.o', 'driver.o', '-o', 'eqquasi_roundtrip'] + lib, ob)
    r, _ = run([os.path.join(ob, 'eqquasi_roundtrip')], qd)
    qout = os.path.join(qd, 'out', 'fault.r.nc')
    with netCDF4.Dataset(qout) as ds:
        print('== 4. EQquasi netcdf_read_on_fault_restart -> netcdf_write_on_fault: dims %s, '
              'shear_strike dtype %s, sha256 %s' % ({k: len(v) for k, v in ds.dimensions.items()},
                                                   ds['shear_strike'].dtype, sha(qout)))

    # 5. EQdyna driver (this checkout's netcdf_io.f90, item7 harness).
    dd = os.path.join(wd, 'eqdyna_stage')
    os.makedirs(dd)
    shutil.copy(qout, os.path.join(dd, 'fault.r.nc'))
    write_nodes(os.path.join(dd, 'nodes.txt'), with_dz=True)
    eb = os.path.join(wd, 'eqdyna_build')
    os.makedirs(eb)
    item7._DRIVER_SRC = EQD_DRIVER
    binary = item7.build_driver(eb)
    run([binary], dd)
    out = np.loadtxt(os.path.join(dd, 'fric_out.txt'), ndmin=2)

    # 6. compare.
    print('== 6. EQdyna netcdf_read_on_fault_eqdyna_restart vs frt.txt* (original float64)')
    print('   %-5s %7s %14s %14s %12s %10s %10s' % ('fault', 'nodes', 'max|ref|', 'max|got-ref|',
                                                  'max rel', 'nonzero', '!=f32(ref)'))
    ok = True
    for ift in range(ntotft):
        ref = np.array(list(nodes[ift].values()))
        got = out[out[:, 0] == ift + 1][:, 2:]
        ref32 = ref.astype(np.float32).astype(np.float64)
        diff = np.abs(got - ref)
        rel = diff / np.maximum(np.abs(ref), 1e-300)
        bad = int(np.count_nonzero(got != ref32))
        print('   %-5d %7d %14.6e %14.6e %12.3e %10d %10d' % (
            ift + 1, len(ref), np.abs(ref).max(), diff.max(), rel[ref != 0].max() if (ref != 0).any() else 0,
            int(np.count_nonzero(ref)), bad))
        for k, name in enumerate(VARS):
            print('      %-17s max|ref| %13.6e  max|got-ref| %13.6e  !=f32 %d' % (
                name, np.abs(ref[:, k]).max(), diff[:, k].max(), int(np.count_nonzero(got[:, k] != ref32[:, k]))))
        ok = ok and bad == 0 and got.shape == ref.shape
    print('float32 unit roundoff 2^-24 = %.3e' % 2.0**-24)
    print('RESULT', 'PASS' if ok else 'FAIL',
          '-- every restart value EQdyna reads back equals float32(original) exactly' if ok else '')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
