# TPV30 on-fault station time series — EQdyna today, dx = 500 m, Fortran

Figure input for `scripts/figures/make_fig2_tpv30_station_timeseries.py`.
Not a test reference, not compared by any gate, not read by `testsys/`.

## Why these files exist at all

`test.reference.results/test.tpv30/` carries only a final-state
`frt.canonical.txt` and a time-independent `fault.dyna.r.nc` (dims
`dip`×`strike` only, no time dimension). There is therefore **no committed
time series** for TPV30, and figure 2 would not be regenerable from the
repository. These 13 files close that gap.

## What run produced them

A fresh `test.tpv30` run of the committed compset, reproducing the frozen
reference exactly:

    export EQDYNAROOT=<checkout>
    create.newcase <rundir> test.tpv30
    cd <rundir> && python3 ./case.setup && mpirun -np 4 eqdyna
    python3 -m testsys.frt_canonical <rundir>

| field | value |
|---|---|
| date | 2026-09-22 |
| host | this box, 4 MPI ranks (`par.nx,ny,nz = 2,1,2`) |
| binary | `bin/eqdyna`, built from the checkout's own `src/fortran` |
| `par.dx` | 500 m (fault grid 81 × 41 = 3321 nodes) |
| `par.dt` | 0.0416667 s, 480 steps, `par.term` = 20 s |
| `par.tpv` | 30 |

**Verified identical to the frozen reference**: this run's
`frt.canonical.txt` vs
`test.reference.results/test.tpv30/frt.canonical.txt` —
`max |diff| = 0.000000e+00` over all 3321 × 22 values. The station series
below are therefore the reference run's own time histories, not a different
run's.

`case.setup.log` is that run's setup output, kept for the geometry
validation lines (`nnx=81, nnz=41, dx=500.0 m`, source
`bFault_Rough_Geometry.tpv29.100m.txt`, native spacing 100 m).

## The files

13 of the 24 spec on-fault stations: exactly those whose strike AND down-dip
coordinates are multiples of 500 m, so that each coincides **exactly** with
the same-named station in
`scec_archive/tpv30/eqdyna-v3.1-100m-2015/`. The other 11 spec stations are
not on the 500 m grid; EQdyna snaps them, and snapped stations are excluded
here rather than compared at an offset.

Format: the SCEC on-fault 8-column layout (t, h-slip, h-slip-rate,
h-shear-stress, v-slip, v-slip-rate, v-shear-stress, n-stress); slip in m,
slip rate in m/s, stresses in MPa.

## One defect found in today's writer, reported not patched

The header block these files carry says `# Time series in 11 columns in
format E15.7` and then writes **8** columns. The 2015 archive's own header
says 8 and writes 8. The data is fine; the header line is wrong. That line
is emitted by the Fortran station writer, which is outside this work's
remit — recorded here so the next reader does not chase it.
