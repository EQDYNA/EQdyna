"""checkInputConsistency.py <- src/fortran/checkInputConsistency.f90.

Refuses EXACTLY what the Fortran refuses: same condition, same message
text, same numbered exit code (src/fortran/errorCodes.f90), at the same
point in the run -- eqdyna3d.f90:104 calls checkInputConsistency right
after read_fault_rough_geometry and before calcLocalShapeFunc/meshgen (i.e.
after every bFile is read, before any mesh or solver work). eqdyna3d.py's
build_solver_state calls check() immediately after
readInputFiles.build_params returns, before its OWN scope guards
(ntotft/friclaw/npx,npy,npz/C_degen) -- none of those three checks below
need anything build_params does not already return, so this is the
earliest point on the Python side too, and it runs for BOTH run_case and
run_case_mpi since both funnel through build_solver_state.

Three checks, ported verbatim (checkInputConsistency.f90:7-19):

  1. C_elastic==0 and C_Q==1              -> ERR_CFG_Q_NEEDS_ELASTIC   (11)
  2. C_Q==1 and rat>1                      -> ERR_CFG_Q_NEEDS_UNIFORM  (12)
  3. output_plastic==1 and C_elastic!=0    -> ERR_CFG_PLASTIC_OUTPUT   (13)

C_Q is not a case-input field on EITHER side: globalvar.f90:104 declares it
`= 0` and no reader (readInputFiles.f90, scripts/case.setup) ever assigns
it from a file -- confirmed by grep, and consistent with
docs/fortran_python_correspondence.md's note on calcQAttenuationCoeff.f90
("C_Q is hardcoded 0 ... and read from no input file"). This module
hardcodes C_Q=0 for the identical reason (module-level `C_Q` below), which
makes checks 1 and 2 genuinely UNREACHABLE through any real case on BOTH
sides. They are still ported, per rule 23 ("the port refuses exactly what
Fortran refuses" -- not just the reachable subset of it), and `check()`
accepts C_Q as an explicit keyword so a future case format that DOES
surface it can exercise both checks without editing this function; see
testsys/regression/test_check_input_consistency.py for how the two
unreachable checks are pinned (direct unit calls, not a case directory) vs.
the one reachable check (output_plastic/C_elastic, real bGlobal.txt on
both binaries).
"""


class InputConsistencyError(RuntimeError):
    """Raised by check(). `.code` is the matching src/fortran/errorCodes.f90
    ERR_CFG_* value, so `python3 -m eqdyna` (see eqdyna3d.main's `_abort`)
    can exit with the SAME number a Fortran run of the same bad config
    exits with (rule 23)."""

    def __init__(self, code, message):
        super().__init__(message)
        self.code = code


# --- 11-19 configuration and parameter consistency (errorCodes.f90:63-65) -
ERR_CFG_Q_NEEDS_ELASTIC = 11   # C_Q=1 requires C_elastic=1
ERR_CFG_Q_NEEDS_UNIFORM = 12   # C_Q=1 requires rat=1.0 (uniform elements)
ERR_CFG_PLASTIC_OUTPUT = 13    # output_plastic=1 requires C_elastic=0
ERR_CFG_NSTRESS_SIGN_INVALID = 15  # bGlobal.txt's station n-stress sign is neither +1 nor -1 (raised by readInputFiles.read_bglobal)
ERR_GEOM_MULTIFAULT_Y_BAD = 32  # a fault's y-plane is not vertical/planar, not a multiple of dy, outside the uniform-y belt, or coincides with another fault's
ERR_GEOM_MULTIFAULT_XZ_BAD = 33  # a fault's x/z extent differs from fault 1's (the shared uniform x/z belt only covers fault 1's box)

# C_Q is never a case-input field (see module docstring) -- hardcoded here
# exactly as globalvar.f90:104 hardcodes its default, never overridden by
# any reader on either side.
C_Q = 0


def check(C_elastic, output_plastic, rat, C_Q=C_Q):
    """checkInputConsistency.f90:7-19, verbatim, same order. Raises
    InputConsistencyError on the first failing check with the Fortran's
    message text unchanged."""
    if C_elastic == 0 and C_Q == 1:
        raise InputConsistencyError(
            ERR_CFG_Q_NEEDS_ELASTIC,
            'Q model (C_Q=1) can only work with the elastic code (C_elastic=1). '
            'Set C_Q=0 or C_elastic=1.')
    if C_Q == 1 and rat > 1:
        raise InputConsistencyError(
            ERR_CFG_Q_NEEDS_UNIFORM,
            'Q model (C_Q=1) can only work with uniform element size; rat must be 1.0.')
    if output_plastic == 1 and C_elastic != 0:
        print(' Now, C_elastic = ', C_elastic, flush=True)   # checkInputConsistency.f90:16
        raise InputConsistencyError(
            ERR_CFG_PLASTIC_OUTPUT,
            'Plastic strains are only output for C_elastic=0. Set output_plastic=0 or C_elastic=0.')


def check_multifault(faults, dy, dis4uniF, dis4uniB, C_degen, tol=1.0e-5):
    """checkInputConsistency.f90's row-17 multi-fault guards, verbatim
    (only the `C_degen == 0.0` branch -- C_degen>3's dipping/wedge-
    degeneration mechanism is a different, pre-existing, single-fault-only
    path with a legitimate fymin != fymax, orthogonal to this work).
    No-op at ntotft==1 (every loop below is over a single fault, and the
    i<j distinctness/xz-match loop does not execute for ntotft<2) --
    bit-identical refusal behaviour to before this function existed.

    `faults`: read_bfaultgeometry's return (list of ntotft dicts with
    fxmin/fxmax/fymin/fymax/fzmin/fzmax)."""
    if C_degen != 0.0:
        return
    f0 = faults[0]
    for i, f in enumerate(faults, start=1):
        if abs(f['fymax'] - f['fymin']) > tol:
            raise InputConsistencyError(
                ERR_GEOM_MULTIFAULT_Y_BAD,
                "checkInputConsistency: fault %d has fymin /= fymax -- only a planar, "
                "vertical fault (single y-plane) is supported; non-planar multi-fault "
                "geometry is out of scope." % i)
        y_over_dy = f['fymin'] / dy
        nearest_int = round(y_over_dy)
        if abs(y_over_dy - nearest_int) > 1.0e-6:
            raise InputConsistencyError(
                ERR_GEOM_MULTIFAULT_Y_BAD,
                "checkInputConsistency: fault %d's y = %.3f is not an integer multiple of "
                "dy = %.3f -- it would fall between mesh node lines and mesh with zero "
                "fault nodes." % (i, f['fymin'], dy))
        if f['fymin'] < -float(dis4uniF) * dy - tol or f['fymin'] > float(dis4uniB) * dy + tol:
            raise InputConsistencyError(
                ERR_GEOM_MULTIFAULT_Y_BAD,
                "checkInputConsistency: fault %d's y = %.3f lies outside the uniform-y mesh "
                "belt [-dis4uniF*dy, +dis4uniB*dy] = [%.3f, %.3f]; widen "
                "par.nuni_y_minus/par.nuni_y_plus or move the fault."
                % (i, f['fymin'], -float(dis4uniF) * dy, float(dis4uniB) * dy))
        if i > 1:
            if (abs(f['fxmin'] - f0['fxmin']) > tol or abs(f['fxmax'] - f0['fxmax']) > tol or
                    abs(f['fzmin'] - f0['fzmin']) > tol or abs(f['fzmax'] - f0['fzmax']) > tol):
                raise InputConsistencyError(
                    ERR_GEOM_MULTIFAULT_XZ_BAD,
                    "checkInputConsistency: fault %d has a different x/z extent than fault 1. "
                    "The shared uniform x/z mesh belt is built from fault 1 box alone, so every "
                    "fault must share it (two parallel faults, same strike extent) -- "
                    "independent per-fault x/z extents are out of scope." % i)
    for i in range(len(faults)):
        for j in range(i + 1, len(faults)):
            if abs(faults[i]['fymin'] - faults[j]['fymin']) < tol:
                raise InputConsistencyError(
                    ERR_GEOM_MULTIFAULT_Y_BAD,
                    "checkInputConsistency: faults %d and %d are both at y = %.3f -- two "
                    "faults must occupy distinct y-planes." % (i + 1, j + 1, faults[i]['fymin']))
