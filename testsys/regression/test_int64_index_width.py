#! /usr/bin/env python3
"""
Regression guard (pathway item 143): the Fortran solver's per-rank node /
element / equation COUNTS, the node and equation IDS, and the OFFSETS into
eqNumIndexArr and stressArr are 64-bit, and the grid->id arithmetic that
feeds them does not wrap past 2^31-1.

WHY THIS EXISTS. Until item 143 `globalvar.f90` declared totalNumOfNodes,
totalNumOfElements, totalNumOfEquations, sizeOfEqNumIndexArr and
sizeOfStressDofIndexArr as integer(kind=4), and meshgen.f90 built the
master-node id as the inline product `nx*ny*nz + ...` with kind=4 operands.
gfortran does not trap integer overflow, so a product past 2,147,483,647
WRAPS to a valid-looking negative or small positive id -- the silent
corruption class the jax port refuses loudly (`backend.check_index_width`).
The offsets wrap first: eqNumStartIndexLoc reaches 12*N at PML nodes and
stressCompIndexArr 21*E at PML elements, and allocInit sizes stressArr as
5*sizeOfEqNumIndexArr -- all past 2^31 at ~1.4e8-1.8e8 elements on ONE
rank, which a fat-memory node can hold; the raw counts follow at 2^31.

WHAT IT CHECKS, without building any mesh (a 2^31-node mesh is ~1 TB):

  1. STATIC: the named declarations in globalvar.f90 are integer(kind=8),
     and meshgen.f90 / assembleGlobalMass.f90 route their id arithmetic
     through the module procedures gridNodeCount / regularNodeId (so the
     inline kind=4 product cannot quietly come back).
  2. MODULE-VARIABLE PATH (driver A): compiled against the REAL
     globalvar.f90, stores exact 64-bit values into the production count
     variables and prints them back; every value must round-trip past 2^31.
     On the pre-143 declarations this COMPILES and prints the wrapped
     values -- the symptom, with numbers.
  3. HELPER PATH (driver B): calls gridNodeCount and regularNodeId with
     (nx, ny, nz) whose product exceeds 2^31-1 and must get the exact
     Python-integer answer, including the LAST node id == nx*ny*nz.
     Pre-143 code has no such procedures: compile error, reported as such.
  4. SELF-TEST (this guard can go red): driver A is compiled against a copy
     of globalvar.f90 with the five count declarations mutated back to
     kind=4; the run MUST reproduce the wrapped values this test exists to
     catch (2400000007 -> -1894967289 etc.). A guard whose red side has
     never been seen tests nothing (papercuts.md).

Cheap (rule 9): three small gfortran builds (globalvar.f90 + one driver),
under 3 s. Fails loudly (rule 2) if gfortran is unavailable -- never
silently skips. No MPI, no case, no solver run.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FSRC = os.path.join(ROOT, 'src', 'fortran')
I32_MAX = 2**31 - 1

WIDE_SCALARS = ('totalNumOfNodes', 'totalNumOfElements', 'totalNumOfEquations',
                'sizeOfEqNumIndexArr', 'sizeOfStressDofIndexArr')
WIDE_ARRAYS = ('eqNumIndexArr', 'stressCompIndexArr', 'eqNumStartIndexLoc',
               'surfaceNodeIdArr', 'nodeElemIdRelation', 'idhist',
               'OffFaultStNodeIdIndex', 'nsmp')

# (nx, ny, nz, nft): the first wraps the NODE count itself (2.4e9 nodes);
# the second keeps N under 2^31 (2.16e8) but wraps the 12*N / 21*E OFFSETS.
GRIDS = ((2000, 2000, 600, 7), (600, 600, 600, 3))

# Driver A: the production module variables, written and read back.
DRIVER_A = r"""
program int64_driver_a
    use globalvar
    implicit none
    integer (kind = 4) :: nx, ny, nz, nft
    integer (kind = 8) :: nodes, elems
    read(*,*) nx, ny, nz, nft
    nodes = int(nx,8)*int(ny,8)*int(nz,8)
    elems = (int(nx,8)-1_8)*(int(ny,8)-1_8)*(int(nz,8)-1_8)
    totalNumOfNodes         = nodes + int(nft,8)
    totalNumOfElements      = elems
    totalNumOfEquations     = 3_8*totalNumOfNodes
    sizeOfEqNumIndexArr     = 12_8*totalNumOfNodes
    sizeOfStressDofIndexArr = 21_8*totalNumOfElements
    write(*,'(a,i0)') 'totalNumOfNodes=',         totalNumOfNodes
    write(*,'(a,i0)') 'totalNumOfElements=',      totalNumOfElements
    write(*,'(a,i0)') 'totalNumOfEquations=',     totalNumOfEquations
    write(*,'(a,i0)') 'sizeOfEqNumIndexArr=',     sizeOfEqNumIndexArr
    write(*,'(a,i0)') 'stressArrSize=',           5*sizeOfEqNumIndexArr
    write(*,'(a,i0)') 'sizeOfStressDofIndexArr=', sizeOfStressDofIndexArr
end program int64_driver_a
"""

# Driver B: the module procedures meshgen/MPI4NodalQuant now route through.
DRIVER_B = r"""
program int64_driver_b
    use globalvar
    implicit none
    integer (kind = 4) :: nx, ny, nz, nft, ix, iy, iz
    read(*,*) nx, ny, nz, nft
    write(*,'(a,i0)') 'gridNodeCount=',  gridNodeCount(nx, ny, nz)
    write(*,'(a,i0)') 'msnode=',         gridNodeCount(nx, ny, nz) + nft
    write(*,'(a,i0)') 'lastNodeId=',     regularNodeId(nx, ny, nz, ny, nz)
    ix = nx/2 + 1; iy = ny/3 + 1; iz = nz/4 + 1
    write(*,'(a,i0)') 'midNodeId=',      regularNodeId(ix, iy, iz, ny, nz)
    write(*,'(a,i0)') 'firstNodeId=',    regularNodeId(1, 1, 1, ny, nz)
end program int64_driver_b
"""


def expected_a(nx, ny, nz, nft):
    nodes = nx*ny*nz
    elems = (nx-1)*(ny-1)*(nz-1)
    return {
        'totalNumOfNodes': nodes + nft,
        'totalNumOfElements': elems,
        'totalNumOfEquations': 3*(nodes + nft),
        'sizeOfEqNumIndexArr': 12*(nodes + nft),
        'stressArrSize': 5*12*(nodes + nft),
        'sizeOfStressDofIndexArr': 21*elems,
    }


def expected_b(nx, ny, nz, nft):
    nodes = nx*ny*nz
    ix, iy, iz = nx//2 + 1, ny//3 + 1, nz//4 + 1
    return {
        'gridNodeCount': nodes,
        'msnode': nodes + nft,
        'lastNodeId': nodes,
        'midNodeId': (ix-1)*ny*nz + (iz-1)*ny + iy,
        'firstNodeId': 1,
    }


def wrap32(v):
    """What a kind=4 store of v prints: two's-complement wrap to 32 bits."""
    v &= 0xFFFFFFFF
    return v - 2**32 if v >= 2**31 else v


def gfortran(args):
    return subprocess.run(['gfortran', '-O0', '-ffree-line-length-none'] + args,
                          capture_output=True, text=True)


def build(tmp, globalvar_src, driver_src, tag):
    """Compile globalvar.f90 (real or mutated) + a driver with plain gfortran.
    Returns (binary_path, None) or (None, compiler_output)."""
    if shutil.which('gfortran') is None:
        raise RuntimeError('gfortran not found on PATH -- this guard cannot run (rule 2: no silent skip)')
    mod_dir = os.path.join(tmp, tag)
    os.makedirs(mod_dir, exist_ok=True)
    gv_obj = os.path.join(mod_dir, 'globalvar.o')
    r = gfortran(['-J', mod_dir, '-c', globalvar_src, '-o', gv_obj])
    if r.returncode != 0:
        return None, r.stderr
    drv = os.path.join(mod_dir, 'driver.f90')
    with open(drv, 'w') as f:
        f.write(driver_src)
    binary = os.path.join(mod_dir, 'driver')
    r = gfortran(['-I', mod_dir, '-J', mod_dir, drv, gv_obj, '-o', binary])
    if r.returncode != 0:
        return None, r.stderr
    return binary, None


def run(binary, grid):
    r = subprocess.run([binary], input=' '.join(str(g) for g in grid) + '\n',
                       capture_output=True, text=True, timeout=30)
    if r.returncode != 0:
        raise RuntimeError(f'driver exited {r.returncode}: {r.stderr}')
    out = {}
    for line in r.stdout.splitlines():
        k, _, v = line.partition('=')
        out[k.strip()] = int(v)
    return out


def check_static(failures):
    gv = open(os.path.join(FSRC, 'globalvar.f90')).read()
    n_checked = 0
    for name in WIDE_SCALARS:
        m = re.search(r'^\s*integer\s*\(\s*kind\s*=\s*(\d)\s*\)\s*::[^!\n]*\b' + name + r'\b', gv, re.M)
        n_checked += 1
        if m is None:
            failures.append(f'static: globalvar.f90 has no scalar declaration for {name}')
        elif m.group(1) != '8':
            failures.append(f'static: globalvar.f90 declares {name} integer(kind={m.group(1)}); must be kind=8')
    joined = re.sub(r'&\s*\n\s*(&)?', ' ', gv)   # join continuation lines
    for name in WIDE_ARRAYS:
        m = re.search(r'^\s*integer\s*\(\s*kind\s*=\s*(\d)\s*\)\s*,\s*allocatable[^!\n]*\b' + name + r'\b', joined, re.M)
        n_checked += 1
        if m is None:
            failures.append(f'static: globalvar.f90 has no allocatable declaration for {name}')
        elif m.group(1) != '8':
            failures.append(f'static: globalvar.f90 declares {name} integer(kind={m.group(1)}); must be kind=8')
    code_only = lambda s: '\n'.join(l.split('!')[0] for l in s.splitlines())
    mg = code_only(open(os.path.join(FSRC, 'meshgen.f90')).read())
    am = code_only(open(os.path.join(FSRC, 'assembleGlobalMass.f90')).read())
    n_gnc = mg.count('gridNodeCount(')
    n_rni = am.count('regularNodeId(')
    if n_gnc < 3:
        failures.append(f'static: meshgen.f90 calls gridNodeCount( {n_gnc} times; expected >= 3 '
                        '(msnode init, createMasterNode, meshGenError)')
    if n_rni != 6:
        failures.append(f'static: assembleGlobalMass.f90 calls regularNodeId( {n_rni} times; expected 6 '
                        '(3 faces x fetch/add)')
    inline = re.findall(r'numxyz\(2\)\s*\*\s*numxyz\(3\)\s*\+', am)
    if inline:
        failures.append(f'static: assembleGlobalMass.f90 still has {len(inline)} inline kind=4 node-id product(s)')
    print(f'  static: {n_checked} declarations checked; gridNodeCount x{n_gnc} in meshgen.f90, '
          f'regularNodeId x{n_rni} in assembleGlobalMass.f90')


def mutate_to_kind4(src_path, dst_path):
    """Revert the five count scalars to kind=4 -- the pre-item-143 code.
    Returns how many declaration statements were changed (0 = nothing to
    mutate, i.e. the source is ALREADY the pre-143 kind=4 code)."""
    s = open(src_path).read()
    n = 0
    for name in WIDE_SCALARS:
        pat = re.compile(r'^(\s*)integer\s*\(\s*kind\s*=\s*8\s*\)(\s*::[^!\n]*\b' + name + r'\b)', re.M)
        s, k = pat.subn(r'\1integer (kind = 4)\2', s)
        n += k
    with open(dst_path, 'w') as f:
        f.write(s)
    return n


def main():
    failures = []
    print('test_int64_index_width: item 143 -- 64-bit counts/ids/offsets past 2^31')
    check_static(failures)
    real_gv = os.path.join(FSRC, 'globalvar.f90')

    with tempfile.TemporaryDirectory(prefix='eqdyna_int64_') as tmp:
        # --- driver A against the real module: the production count variables.
        bin_a, err = build(tmp, real_gv, DRIVER_A, 'real_a')
        if bin_a is None:
            failures.append('driver A (module variables) failed to compile against the real globalvar.f90:\n'
                            + err[-1500:])
        else:
            n_vals = 0
            for grid in GRIDS:
                exp = expected_a(*grid)
                got = run(bin_a, grid)
                for k, v in exp.items():
                    n_vals += 1
                    if got.get(k) != v:
                        failures.append(f'module var, grid {grid}: {k} = {got.get(k)}, expected {v} '
                                        f'(kind=4 wrap of the expected value is {wrap32(v)})')
                big = [k for k, v in exp.items() if v > I32_MAX]
                print(f'  module vars, grid {grid}: {len(exp)} values, {len(big)} past 2^31-1 ({", ".join(big)})')
            print(f'  module vars: {n_vals} values compared against exact Python integers')

        # --- driver B against the real module: the id-arithmetic helpers.
        bin_b, err = build(tmp, real_gv, DRIVER_B, 'real_b')
        if bin_b is None:
            failures.append('driver B (gridNodeCount/regularNodeId) failed to compile against the real '
                            'globalvar.f90 -- the 64-bit id helpers are missing:\n' + err[-800:])
        else:
            n_vals = 0
            for grid in GRIDS:
                exp = expected_b(*grid)
                got = run(bin_b, grid)
                for k, v in exp.items():
                    n_vals += 1
                    if got.get(k) != v:
                        failures.append(f'helpers, grid {grid}: {k} = {got.get(k)}, expected {v} '
                                        f'(kind=4 wrap of the expected value is {wrap32(v)})')
                big = [k for k, v in exp.items() if v > I32_MAX]
                print(f'  helpers, grid {grid}: {len(exp)} values, {len(big)} past 2^31-1 ({", ".join(big)})')
            print(f'  helpers: {n_vals} values compared against exact Python integers')

        # --- self-test: the pre-143 declarations MUST wrap, visibly.
        mutated = os.path.join(tmp, 'globalvar.k4.f90')
        n_mut = mutate_to_kind4(real_gv, mutated)
        if n_mut == 0:
            failures.append('self-test: no kind=8 count declaration to mutate -- globalvar.f90 is already '
                            'the pre-143 kind=4 code (see the static and module-variable findings above)')
        else:
            bin_4, err4 = build(tmp, mutated, DRIVER_A, 'k4')
            if bin_4 is None:
                failures.append('self-test: mutated kind=4 globalvar failed to compile:\n' + err4[-1500:])
            else:
                grid = GRIDS[0]
                exp = expected_a(*grid)
                got4 = run(bin_4, grid)
                wrapped = {k: got4[k] for k in WIDE_SCALARS if got4[k] != exp[k]}
                if len(wrapped) != len(WIDE_SCALARS):
                    failures.append(f'self-test: with {n_mut} declaration statements reverted to kind=4 only '
                                    f'{len(wrapped)} of {len(WIDE_SCALARS)} counts wrapped -- the guard cannot go red')
                for k, v in wrapped.items():
                    if v != wrap32(exp[k]):
                        failures.append(f"self-test: {k} wrapped to {v}, not the two's-complement {wrap32(exp[k])}")
                print(f'  self-test: {n_mut} declaration statements mutated to kind=4; '
                      f'{len(wrapped)}/{len(WIDE_SCALARS)} counts wrapped, e.g. totalNumOfNodes '
                      f'{exp["totalNumOfNodes"]} -> {got4["totalNumOfNodes"]}')

    if failures:
        print('FAIL test_int64_index_width:')
        for f in failures:
            print('  - ' + f)
        return 1
    print('SUCCESS test_int64_index_width')
    return 0


if __name__ == '__main__':
    sys.exit(main())
