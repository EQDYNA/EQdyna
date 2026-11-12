"""Set dx, term and partition in a scratch copy of a tpv36/37 case. usage: set_case.py <case_dir> <dx> <term_s> <nx> <ny> <nz>
Edits only the scratch case directory's tpv36_37_common.py; case_input/ is never touched."""
import re, sys
d, dx, term, nx, ny, nz = sys.argv[1:7]
p = d + "/tpv36_37_common.py"
s = open(p).read()
subs = [(r"(?m)^(\s*par\.dx\s*=\s*)[0-9.]+", r"\g<1>%s." % dx),
        (r"(?m)^(\s*par\.term\s*=\s*)[0-9.]+", r"\g<1>%s" % term),
        (r"(?m)^(\s*par\.nx\s*=\s*)[0-9]+", r"\g<1>%s" % nx),
        (r"(?m)^(\s*par\.ny\s*=\s*)[0-9]+", r"\g<1>%s" % ny),
        (r"(?m)^(\s*par\.nz\s*=\s*)[0-9]+", r"\g<1>%s" % nz)]
for pat, rep in subs:
    s, n = re.subn(pat, rep, s)
    assert n == 1, (pat, n)
open(p, "w").write(s)
for line in s.splitlines():
    if re.match(r"\s*par\.(dx|term|nx|ny|nz|dt|C_degen|insertFaultType|dip)\s*=", line):
        print(line.strip())
