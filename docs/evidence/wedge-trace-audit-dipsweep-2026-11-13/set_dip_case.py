"""Scratch copy of test.tpv36 with ONLY the dip changed (plus what the dip forces, see README).
usage: set_dip_case.py <case_dir> <dip_deg_int> <dx> <term_s> <nx> <ny> <nz> [dt_factor=0.5]
Edits tpv36_37_common.py in the scratch case dir only (case_input/ is never touched):
  par.dip            -> dip
  par.zmin           -> -(faultWidth*sin(dip) + 7.753 km): model bottom held 7.753 km below the fault bottom
                        (equals the shipped -15 km exactly at 15 deg; the 28 km fault does not fit above -15 km at >=45 deg)
  par.dt             -> dt_factor*min(dy,dz)/vp (default 0.5; identical to the shipped 0.5*dz/vp at 15 deg, where dz < dy; 45 deg is unstable at 0.5 and uses 0.25)
  par.dx / term / nx / ny / nz."""
import os, re, sys
d, dip, dx, term, nx, ny, nz = sys.argv[1:8]
dtf = sys.argv[8] if len(sys.argv) > 8 else '0.5'
p = d + "/tpv36_37_common.py"
s = open(p).read()
subs = [(r"(?m)^(\s*par\.dip\s*=\s*)15\b", r"\g<1>%s" % dip),
        (r"(?m)^(\s*par\.zmin, par\.zmax = ).*$",
         r"\g<1>-(faultWidth*sin(abs(par.dip)/180.*pi) + (15.0e3 - 28.e3*sin(15./180.*pi))), 0.0e3"),
        (r"(?m)^(\s*par\.dt\s*=\s*)0\.5\*par\.dz/par\.vp", r"\g<1>%s*min(par.dy, par.dz)/par.vp" % dtf),
        (r"(?m)^(\s*par\.dx\s*=\s*)[0-9.]+", r"\g<1>%s." % dx),
        (r"(?m)^(\s*par\.term\s*=\s*)[0-9.]+", r"\g<1>%s" % term),
        (r"(?m)^(\s*par\.nx\s*=\s*)[0-9]+", r"\g<1>%s" % nx),
        (r"(?m)^(\s*par\.ny\s*=\s*)[0-9]+", r"\g<1>%s" % ny),
        (r"(?m)^(\s*par\.nz\s*=\s*)[0-9]+", r"\g<1>%s" % nz)]
if os.environ.get("DDSIGN") == "-1":   # EXTRA probe (not part of the owner's sweep): reverse the sign of the initial along-dip shear -> normal-faulting sense
    subs.append((r"(?m)^(\s*par\.on_fault_vars\[iz,ix,49\]\s*=\s*)0\.00212e6", r"\g<1>-0.00212e6"))
for pat, rep in subs:
    s, n = re.subn(pat, rep, s)
    assert n == 1, (pat, n)
open(p, "w").write(s)
for line in s.splitlines():
    if re.match(r"\s*par\.(dx|term|nx|ny|nz|dt|C_degen|insertFaultType|dip|zmin, par\.zmax)\s*[=,]", line):
        print(line.strip())
