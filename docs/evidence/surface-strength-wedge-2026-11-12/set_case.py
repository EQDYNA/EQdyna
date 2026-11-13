"""Apply ONE surface-strength change plus dx/term/partition to a scratch case dir (never case_input/).
usage: set_case.py <kind> <case_dir> <dx> <nx> <ny> <nz> <term> <wedge_params_src>
kind: w12c4  TPV12 wedge variant (row149 params file), fric_cohesion 0.2e6 -> 4.0e6
      i12c4  TPV12 insertion (case copy), fric_cohesion 0.2e6 -> 4.0e6
      i12c02 TPV12 insertion (case copy), unchanged cohesion (baseline control)
      w12c02 TPV12 wedge variant, unchanged cohesion (decomposition control)
      t36c02 test.tpv36, cohesionSlope 0.0005e6 -> 0.000025e6 (trace 4 -> 0.2 MPa)
      t36c4  test.tpv36, unchanged cohesion (baseline control)"""
import re, sys
kind, d, dx, nx, ny, nz, term, wsrc = sys.argv[1:9]

def sub(s, pat, rep, tag):
    s, n = re.subn(pat, rep, s)
    if n != 1:
        raise SystemExit("%s: pattern %r matched %d times" % (tag, pat, n))
    return s

if kind.startswith(("w12", "i12")):
    src = wsrc if kind.startswith("w12") else d + "/user_defined_params.py"
    s = open(src).read()
    s = sub(s, r"(?m)^par\.dx *= *[0-9.]+", "par.dx = %s." % dx, "dx")
    s = sub(s, r"(?m)^par\.term *= *[0-9.]+", "par.term = %s" % term, "term")
    s = sub(s, r"(?m)^par\.nx, par\.ny, par\.nz *= *[0-9]+, *[0-9]+, *[0-9]+", "par.nx, par.ny, par.nz = %s, %s, %s" % (nx, ny, nz), "part")
    if kind.endswith("u4"):   # UNIFORM 4 MPa: measured degenerate (no nucleation); kept only as the 500 m evidence run
        s = sub(s, r"(?m)^par\.fric_cohesion     = 0\.2e6$", "par.fric_cohesion     = 4.0e6  # SURFACE-STRENGTH TEST: was 0.2e6", "coh")
    elif kind.endswith("c4"):  # surface-only: 4 MPa at the trace, linear to the unchanged 0.2 MPa at 8 km down dip (nucleation at 10.5-13.5 km untouched)
        s = sub(s, r"(?m)^(    par\.on_fault_vars\[iz,ix,4\] *= *)par\.fric_cohesion$",
                r"\g<1>par.fric_cohesion + (4.0e6 - par.fric_cohesion)*(0.0 if downDipDistance >= 8.0e3 else 1.0 - downDipDistance/8.0e3)"
                r"  # SURFACE-STRENGTH TEST: 4 MPa at trace, linear to 0.2 MPa at 8 km down dip", "coh")
    elif kind not in ("i12c02", "w12c02"):
        raise SystemExit("unknown kind " + kind)
    open(d + "/user_defined_params.py", "w").write(s)
    for line in s.splitlines():
        if re.match(r"par\.(dx|dy|term|C_degen|insertFaultType|dt|fric_cohesion) |par\.nx, par\.ny|\s+par\.on_fault_vars\[iz,ix,4\]", line):
            print(line)
elif kind.startswith("t36"):
    p = d + "/tpv36_37_common.py"
    s = open(p).read()
    for k, v in (("dx", dx + "."), ("term", term), ("nx", nx), ("ny", ny), ("nz", nz)):
        s = sub(s, r"(?m)^(\s*par\.%s\s*=\s*)[0-9.]+" % k, r"\g<1>%s" % v, k)
    open(p, "w").write(s)
    u = open(d + "/user_defined_params.py").read()
    if kind == "t36c02":
        u = sub(u, r"cohesionSlope=0\.0005e6\)", "cohesionSlope=0.000025e6)  # SURFACE-STRENGTH TEST: was 0.0005e6", "coh")
    elif kind != "t36c4":
        raise SystemExit("unknown kind " + kind)
    open(d + "/user_defined_params.py", "w").write(u)
    for line in (s + u).splitlines():
        if re.match(r"\s*par\.(dx|term|nx|ny|nz|dt|C_degen) *=|.*cohesionSlope=", line):
            print(line.strip())
else:
    raise SystemExit("unknown kind " + kind)
