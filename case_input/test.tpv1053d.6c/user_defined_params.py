#! /usr/bin/env python3

from defaultParameters import parameters

par = parameters()
par.term = 5.
par.nx, par.ny, par.nz = 1,2,3

print(dir(par))

