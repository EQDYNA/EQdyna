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
