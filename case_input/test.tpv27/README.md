# test.tpv27

SCEC TPV27: the test.tpv26 benchmark (planar strike-slip fault, same
geometry, same nucleation, same stations) with off-fault Drucker-Prager
viscoplasticity added. Spec: TPV26_27_Description_v13, Part 5 p.5:
"The material properties are the only difference between the elastic
benchmark (TPV26), and the viscoplastic benchmark (TPV27)." See
`case_input/test.tpv26/README.md` for everything that is unchanged.

## Relationship to test.tpv26

This compset is test.tpv26 plus a "--- TPV27 off-fault plasticity ---" block
in `user_defined_params.py` (`par.C_elastic=0`, Drucker-Prager cohesion/bulk
friction, viscoplastic relaxation time, deviatoric-stress depth taper,
plastic-strain output window) -- the same machinery wired for test.tpv30,
confirmed present by reading `meshgen.f90:setPlasticStress` /
`func_lib.f90:devStrDepthTaper` and the matching python port
(`eqdyna3d.py`/`readInputFiles.py`) before writing this case, not assumed
from the TPV30 name alone.

## Viscoplastic parameters (spec Part 5, implied material table)

| parameter | value | case input |
|---|---|---|
| cohesion c | 1.36 MPa | `par.coheplas` |
| bulk friction nu | 0.1934 | `par.bulk` |
| viscoplastic relaxation time Tv | 0.03 s | `par.viscoplasticRelaxTime` |
| deviatoric-stress depth taper | 15-20 km | `par.devStrTaperDepthStart/End` |
| fluid pressure Pf | hydrostatic | `par.gamar = 0` |

`par.str1ToFaultAngle`/`par.devStrToStrVertRatio` (56.708795263101166 deg /
0.1842011609355381) are derived from this case's own (b11, b33, b13) via
`R*cos(2a) = b11-1`, `R*sin(2a) = -b13` -- the same relation verified to
reproduce test.tpv30's own already-committed values before being applied
here.

## Resolution tier

Same as test.tpv26: gated at `par.dx=500` m; the spec's 100 m/50 m request
recorded (not run) in `testsys/e2e/full_specs.py['test.tpv27']`.

## Stations

Identical station grid to test.tpv26 (spec Part 8/9, shared format).
