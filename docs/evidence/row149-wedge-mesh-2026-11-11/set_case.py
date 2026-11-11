"""Set dx and the partition in a copy of the wedge parameter file. usage: set_case.py <src> <dst> <dx> <nx> <ny> <nz>"""
import re, sys
src, dst, dx, nx, ny, nz = sys.argv[1:7]
s = open(src).read()
s, n1 = re.subn(r"(?m)^par\.dx *= *[0-9.]+", "par.dx = %s." % dx, s)
s, n2 = re.subn(r"(?m)^par\.nx, par\.ny, par\.nz *= *[0-9]+, *[0-9]+, *[0-9]+", "par.nx, par.ny, par.nz = %s, %s, %s" % (nx, ny, nz), s)
assert (n1, n2) == (1, 1), (n1, n2)
open(dst, "w").write(s)
for line in s.splitlines():
    if re.match(r"par\.(dx|dy|term|C_degen|insertFaultType|dt|nuni_y_plus|fymin) |par\.nx, par\.ny", line):
        print(line)
