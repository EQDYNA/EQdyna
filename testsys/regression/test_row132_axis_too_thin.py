#! /usr/bin/env python3
"""
Regression guard (board row 132(1)): a rank whose local 1D slice on x/y/z
holds fewer than 2 nodes is refused LOUDLY, in both languages, instead of
setSurfaceStation's elseif chain silently dropping a station it should
have owned (see docs/notes/NOTES_row131.md for the full derivation: when a
rank's local nx==1, ix==1 is tested before ix==nx in the elseif chain
and always wins, so the ALWAYS-a-candidate high-edge branch becomes
unreachable for that rank).

WHY A STANDALONE PROBE, NOT A FULL CASE RUN: every real case in this repo
(tpv8, tpv10, tpv36, tpv37, drv.a6) has a mesh floor -- even coarsened to
the largest dx/nPML this suite own case.setup will accept, tpv8 own
global node counts bottom out at 17 (x), 25 (y), 9 (z) (measured directly,
docs/notes/NOTES_row131.md) -- all well above the shared box 8-rank cap,
so no REAL case run driven through case.setup + mpirun can put fewer
than 2 nodes on any rank axis within that budget. getLocalOneDimCoorArrAndSize
(src/fortran/meshgen.f90) and partition_1d/check_partition_1d
(src/python/eqdyna/meshgen.py) are pure functions of
(fault span, dx, xmin/xmax, rat, nPML, npx) and (global_size, num_mpi)
respectively -- this probe calls the REAL production routines directly
with a small, hand-picked, but real and internally-consistent set of those
inputs, rather than mocking or reimplementing the arithmetic.

THE DECOMPOSITION UNDER TEST: dx=1, fault span (fltxyz(2,1,1)-fltxyz(1,1,1))=1,
xmin=-2, xmax=2, rat=1.025, nPML=1, npx=8. Hand-derived and cross-checked
against Python own partition_1d (verbatim port, rule 23): the global 1D
line has exactly 7 nodes, split across 8 ranks with numOfNodesPerMPI=1,
residualNumOfNodes=6 -- ranks 0 AND 1 both get a 1-node local slice.
Rank 0 is the global low edge; rank 1 is NOT an edge (ranks 0 and 7 are the
two edges of an 8-rank line) -- exactly the MIDDLE-rank shape this row
names, where setSurfaceStation elseif chain silently breaks. Verified
by hand (docs/notes/NOTES_row131.md) and mechanically here via Python
check_partition_1d(7, 8), which raises citing rank 1 holds 1 node(s)
explicitly (not just rank 0).

FORTRAN: writes a tiny standalone program to a temp dir, sets the module
globals above directly (use globalvar), and calls
getLocalOneDimCoorArrAndSize with MPIXyzId=1 (rank 1, the middle rank)
and dimId=1 (x axis) -- the exact call meshgen own node loop makes,
just invoked directly instead of through the full mesh-generation driver.
No mpirun needed: MPI_Initialized is false before MPI_Init is ever
called, so abortRun (src/fortran/errorCodes.f90) falls through to a plain
call exit(code). Linked against the SAME globalvar.o/errorCodes.o/
meshgen.o a normal make produces (build_eqdyna, shared with
test_row127/test_row120) -- never a reimplementation.

PYTHON: calls meshgen.check_partition_1d(7, 8) directly -- the identical
(global_size, num_mpi) pair -- and asserts it raises, citing rank 1 (not
only rank 0) by name.

Cost: reuses the incremental Fortran build test_row127/test_row120 already
pay for; one extra ~20-line program compile+link (no mpirun); one Python
call. Sub-second beyond the shared build.
"""
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FSRC = os.path.join(ROOT, "src", "fortran")
PYSRC = os.path.join(ROOT, "src", "python")
sys.path.insert(0, PYSRC)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_row127_station_ownership as row127  # reuses build_eqdyna, no reimplementation

ERR_MPI_AXIS_TOO_THIN = 53  # src/fortran/errorCodes.f90 -- kept in sync by test_stop_exit_status.py

PROBE_SOURCE = """\
! Standalone probe for row 132(1) -- see this file module docstring for
! the full derivation of every literal constant below.
program probeAxisTooThin
    use globalvar
    implicit none
    integer (kind = 4) :: globalOneDimCoorArrSize, numOfNodesWithUniformGridsize
    integer (kind = 4) :: frontEdgeNodeId, localOneDimCoorArrSize, globalOneDimCoorArrOutSize
    real (kind = dp) :: localOneDimCoorArr(10000), modelBoundCoor(3,2), globalOneDimCoorArrOut(10000)

    allocate(fltxyz(2,4,1))
    fltxyz(1,1,1) = 0.0d0
    fltxyz(2,1,1) = 1.0d0
    dx   = 1.0d0
    rat  = 1.025d0
    nPML = 1
    xmin = -2.0d0
    xmax = 2.0d0
    npx  = 8

    call getLocalOneDimCoorArrAndSize(globalOneDimCoorArrSize, numOfNodesWithUniformGridsize, &
        frontEdgeNodeId, 1, localOneDimCoorArrSize, localOneDimCoorArr, modelBoundCoor, 1, &
        globalOneDimCoorArrOut, globalOneDimCoorArrOutSize)

    write(*,*) "NO ABORT (bug): globalOneDimCoorArrSize=", globalOneDimCoorArrSize, &
        " localOneDimCoorArrSize=", localOneDimCoorArrSize
end program probeAxisTooThin
"""


def run_fortran_probe():
    """Build+run the standalone probe. Returns (returncode, stdout+stderr)."""
    row127.build_eqdyna()  # ensures FSRC/*.o and *.mod are fresh; shared, incremental
    with tempfile.TemporaryDirectory() as d:
        src = os.path.join(d, "probeAxisTooThin.f90")
        exe = os.path.join(d, "probeAxisTooThin")
        with open(src, "w") as f:
            f.write(PROBE_SOURCE)
        # meshgen.o (one compilation unit with many subroutines) references
        # symbols from other files even though this probe calls only
        # getLocalOneDimCoorArrAndSize -- the linker must resolve every
        # external reference in the .o, not just the one we invoke. Link
        # every object the real `eqdyna` build produces except eqdyna3d.o
        # (its own PROGRAM would collide with this probe's).
        objs = sorted(
            os.path.join(FSRC, n) for n in os.listdir(FSRC)
            if n.endswith('.o') and n != 'eqdyna3d.o')
        missing = [o for o in objs if not os.path.exists(o)]
        if missing:
            raise RuntimeError("probe link: missing %r after build_eqdyna()" % missing)
        r = subprocess.run(
            ["mpif90", "-I", FSRC, "-o", exe, src] + objs +
            ["-L/usr/lib/x86_64-linux-gnu", "-lnetcdf", "-lnetcdff"],
            capture_output=True, text=True, timeout=120)
        if r.returncode != 0:
            raise RuntimeError("probe compile/link failed (exit %d):\n%s"
                               % (r.returncode, (r.stdout + r.stderr)[-3000:]))
        p = subprocess.run([exe], capture_output=True, text=True, timeout=60)
        return p.returncode, p.stdout + p.stderr


def main():
    checks = []
    try:
        rc, out = run_fortran_probe()
    except RuntimeError as e:
        print("FAIL test_row132_axis_too_thin")
        print(" -", e)
        return 1

    checks.append(("fortran probe exits ERR_MPI_AXIS_TOO_THIN (%d)" % ERR_MPI_AXIS_TOO_THIN,
                   rc == ERR_MPI_AXIS_TOO_THIN, True))
    checks.append(("fortran probe prints a FATAL block",
                   "EQdyna: FATAL" in out, True))
    checks.append(("fortran probe names the thin-axis reason",
                   "local 1D slice" in out or "node(s)" in out, True))
    if rc != ERR_MPI_AXIS_TOO_THIN or "EQdyna: FATAL" not in out:
        print("  fortran probe output:")
        print("   ", out.replace("\n", "\n    "))

    from eqdyna import meshgen
    try:
        meshgen.check_partition_1d(7, 8)
        py_raised, py_msg = False, ""
    except ValueError as e:
        py_raised, py_msg = True, str(e)
    checks.append(("python check_partition_1d(7, 8) raises", py_raised, True))
    checks.append(("python refusal names the MIDDLE rank (rank 1), not only rank 0",
                   "rank 1 holds 1 node(s)" in py_msg, True))
    if not py_raised:
        print("  python check_partition_1d(7, 8) did not raise -- port gap")

    for label, got, expect_pass in checks:
        ok = (got is True) == expect_pass
        print(("PASS" if ok else "FAIL") + "  " + label)
    fails = [label for label, got, expect_pass in checks if (got is True) != expect_pass]

    print("%s test_row132_axis_too_thin: %d check(s); %d failure(s)"
          % ("FAIL" if fails else "SUCCESS", len(checks), len(fails)))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
