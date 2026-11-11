"""Derive the wedge-mesh scratch parameter file from test.tpv12's. usage: make_params.py <src> <dst>"""
import sys
src = open(sys.argv[1]).read()
def rep(a, b):
    global src
    assert src.count(a) == 1, (a, src.count(a))
    src = src.replace(a, b)
rep("#! /usr/bin/env python3\n", "#! /usr/bin/env python3\n\n# ROW 149 SCRATCH VARIANT (evidence only, not a registered case): test.tpv12 with the\n# 60-degree fault meshed by WEDGE DEGENERATION (par.C_degen = 60, par.insertFaultType = 0)\n# instead of fault insertion. Lines changed from case_input/test.tpv12 are tagged WEDGE.\n# par.dx below is rewritten per run (500 / 200 / 150).\n")
rep("par.fzmin, par.fzmax = -15.0e3, 0.0e3\n", "par.fzmin, par.fzmax = -15.0e3, 0.0e3\nFAULT_WIDTH = 15.0e3  # WEDGE: down-dip width, captured before mod4dip rescales par.fzmin\npar.fymin, par.fymax = 0.0, FAULT_WIDTH*cos(abs(par.dip)/180.*pi)  # WEDGE: the fault spans grid y\n")
rep("par.dy = par.dx\n", "par.dy = par.dx*cos(abs(par.dip)/180.*pi)  # WEDGE: dz/dy = tan(dip) puts the fault through grid nodes\npar.nuni_y_plus  = round(FAULT_WIDTH/par.dx) + 5  # WEDGE: uniform y belt covers the fault (test.tpv36 rule)\npar.nuni_y_minus = 5\n")
rep("par.term = 8.   # spec: 0-8 s after nucleation\n", "par.term = 4.5\n")
rep("par.insertFaultType = 1  # insert planar dipping fault (test.tpv10's path)\n", "par.C_degen = par.dip    # WEDGE: degenerate hexahedra to wedges along the fault\npar.insertFaultType = 0  # WEDGE: no fault insertion\n")
rep("print(' ')\nprint('DIP: dz is    ', par.dz)", "par.dt = 0.5*min(par.dy, par.dz)/par.vp  # WEDGE: dy = dx/2 is now the smallest edge (insertion: 0.5*dz/vp)\nprint(' ')\nprint('DIP: dz is    ', par.dz)")
open(sys.argv[2], "w").write(src)
