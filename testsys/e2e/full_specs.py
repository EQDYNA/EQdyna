#! /usr/bin/env python3
"""
FULL_SPECS -- SCEC spec-resolution/spec-duration overrides for the e2e
"full" tier (run_e2e_full.py). PROJECT_RULES.md rule 2 (no invented data)
and rule 6 (every performance/spec number carries provenance): every entry
below cites the exact SCEC cvws description page and PDF quote it came
from. Fetched 2026-09-14 into scratch/specs/ (curl + pdftotext -layout);
see that directory for the raw PDFs/text this was read from.

Design contract (per coordinator instruction): the full tier NEVER forks a
*_full compset. It runs the SAME case_input/<case> tree the fast e2e tier
uses; the only difference is a mechanical rewrite of par.dx/par.nx,ny,nz/
par.term in the generated user_defined_params.py after create.newcase,
exactly like testsys/perf/run_scaling.py's make_case() does for rank
counts. There is no test.tpv8_full/ or similar directory anywhere.

nx,ny,nz below is the 16-rank decomposition, reusing
testsys/perf/run_scaling.py's DECOMP[16] = (4, 2, 2) convention (this file
does not invent a separate decomposition scheme).
"""

# case -> {dx (m), term (s), (nx,ny,nz), citation}
FULL_SPECS = {
    'test.tpv29': dict(
        dx=50.0, term=20.0, decomp=(4, 1, 4),
        source=("TPV29_30_Description_v06, Part 8: 50 m preferred / 100 m "
                "acceptable; 0-20 s. NOTE: 50 m is ~119 M elements — an HPC "
                "job (bundle at scratch/tpv29/hpc50m, 1024 ranks ~0.7 h). "
                "ny=1 keeps the fault plane off MPI partitions."),
    ),
    'test.tpv8': dict(
        dx=100., term=15., nx=4, ny=2, nz=2,
        citation=(
            'SCEC TPV8/TPV9 description (scratch/specs/TPV8_desc.pdf via '
            'https://strike.scec.org/cvws/tpv89docs.html -> '
            'download/TPV8_forwebsite.pdf), p.7 "Computations should be run '
            'using the following element-size/node-spacing ... 100m" '
            '(150m only if memory/processor-constrained); p.7 "for the '
            'times 0.0 to 15.0 seconds after nucleation."'
        ),
    ),
    'test.tpv10': dict(
        dx=100., term=15., nx=4, ny=2, nz=2,
        citation=(
            'SCEC TPV10/TPV11 description (scratch/specs/TPV10_11_desc.pdf '
            'via https://strike.scec.org/cvws/tpv10_11docs.html -> '
            'download/TPV10_11_Description_v7.pdf), item 2 "Computations '
            'should be run using 100 m node spacing on the fault plane"; '
            'item 3 "time series results for the times 0.0 to 15.0 seconds '
            'after nucleation."'
        ),
    ),
    'test.tpv104': dict(
        dx=50., term=12., nx=4, ny=2, nz=2,
        citation=(
            'SCEC TPV103/TPV104 description page '
            'https://strike.scec.org/cvws/tpv103_104docs.html: '
            '"The recommended element-size or node-spacing is 50 meters." '
            'Duration from scratch/specs/TPV103_104_desc.pdf ('
            'download/SCEC_validation_slip_law.pdf): "Report the complete '
            'time histories from t = 0 to t = 12 s".'
        ),
    ),
}

# case -> reason it is excluded from the full tier (rule 2: never invent a
# number to fill a gap -- exclude and say so instead).
EXCLUDED = {
    'test.tpv1053d': (
        'SCEC TPV105-3D description+file-formats PDFs (scratch/specs/'
        'TPV105_3D_desc.txt, TPV105_3D_formats.txt, fetched from '
        'https://strike.scec.org/cvws/tpv105_3D_docs.html) specify duration '
        '(0 to 15 s, "Time Histories of Fields...") but NEVER state a '
        'required element-size/node-spacing -- the file-formats PDF only '
        'has a free-text "Node spacing or element size" field for the '
        'submitter to report, not a target value. No number to rewrite '
        'par.dx to without inventing one; excluded until SCEC publishes one '
        'or an owner supplies a citation.'
    ),
    'test.drv.a6': (
        'Internal case, no SCEC spec. No published full-resolution run '
        '(different dx/term than the gated test config) discoverable in '
        'README.md or pastReleaseNotes.md -- both only document the case at '
        'its current test parameters (dx=500m, term=5s). Excluded rather '
        'than invented.'
    ),
    'test.meng2023a': (
        'Internal case, no SCEC spec, no published full-resolution figure '
        'in README.md/pastReleaseNotes.md beyond the current test config '
        '(dx=400m, term=5s). Excluded rather than invented.'
    ),
    'test.meng2023cb': (
        'Same as test.meng2023a (paired case, same provenance gap). '
        'Excluded rather than invented.'
    ),
}
