#! /usr/bin/env python3

from defaultParameters import parameters

par = parameters()
par.nx, par.ny, par.nz = 1,2,3
print(dir(par))

