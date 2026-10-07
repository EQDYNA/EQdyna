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
        # Keys must be nx/ny/nz and citation -- run_e2e_full.py reads
        # spec['nx'], spec['ny'], spec['nz'] and spec['citation'] directly.
        # This entry used decomp=(...)/source=(...) and, being first in the
        # dict, made the whole full tier die with KeyError before any case
        # launched.
        dx=50.0, term=20.0, nx=4, ny=1, nz=4,
        citation=("TPV29_30_Description_v06, Part 8: 50 m preferred / 100 m "
                "acceptable; 0-20 s. NOTE: 50 m is ~119 M elements — an HPC "
                "job (bundle at scratch/tpv29/hpc50m, 1024 ranks ~0.7 h). "
                "ny=1 keeps the fault plane off MPI partitions. Geometry at "
                "this dx comes from the compset's shipped 50 m surface "
                "(bFault_Rough_Geometry.tpv29.50m.txt, an exact decimation of "
                "the official 25 m data); no download step is needed. Before "
                "v5.6.0 only a 100 m surface shipped and this entry could not "
                "be set up at all."),
    ),
    'test.tpv30': dict(
        # rule 17 step 5: recorded, not run. Same fault/geometry/nucleation
        # citation as test.tpv29 (identical surface, identical Part 5/6);
        # the ONLY thing TPV30 adds is Part 7 (Drucker-Prager viscoplasticity,
        # already wired via par.viscoplasticRelaxTime/devStrTaperDepthStart/
        # End -- item 24(b)/(c)). Actually running this needs the 50 m
        # bFault_Rough_Geometry file copied into case_input/test.tpv30/ first
        # (only the 100 m file is shipped there today; see that compset's
        # README) -- a scheduling/storage decision, not a spec gap.
        dx=50.0, term=20.0, nx=4, ny=1, nz=4,
        citation=("TPV29_30_Description_v06, Part 3 (p.15) 'Please submit "
                "results using 50 m node spacing on the fault plane. If you "
                "are unable to run the simulation with 50 m node spacing, "
                "then it is OK to use 100 m node spacing.'; 'Run the model "
                "for times from 0.0 to 20.0 seconds after nucleation.' "
                "ny=1 keeps the fault plane off MPI partitions, same as "
                "test.tpv29's own full-tier entry."),
    ),
    'test.tpv26': dict(
        # rule 17 step 5: recorded, not run. requireFaultGeometryResolution
        # admits only dx=500 m today (the gate tier); no 100 m/50 m geometry
        # file is shipped in case_input/test.tpv26/ -- a scheduling/storage
        # decision, not a spec gap (this case's fault is planar, so the full
        # tier needs only a finer on_fault_vars grid, not a downloaded
        # surface the way tpv29/30's rough fault does).
        dx=100.0, term=13.0, nx=4, ny=1, nz=4,
        citation=("TPV26_27_Description_v13, Part 3 p.9: 'We request that "
                "you run each of these two benchmarks using two resolutions: "
                "100 meter resolution, and 50 meter resolution... If you are "
                "unable to run the simulation with 50 m node spacing, then "
                "it is OK to omit the 50 m case.' 'Run the model for times "
                "from 0.0 to 13.0 seconds after nucleation.' ny=1 keeps the "
                "fault plane off MPI partitions, same convention as every "
                "other planar-fault full-tier entry here."),
    ),
    'test.tpv27': dict(
        # Same geometry/term citation as test.tpv26 (spec Part 5 p.5:
        # "material properties are the only difference" between the two
        # benchmarks) -- the Drucker-Prager machinery is already wired via
        # par.viscoplasticRelaxTime/devStrTaperDepthStart/End (confirmed
        # present from test.tpv30, row 150 investigation).
        dx=100.0, term=13.0, nx=4, ny=1, nz=4,
        citation=("TPV26_27_Description_v13, Part 3 p.9 (TPV26 and TPV27 "
                "share the resolution/term request -- see test.tpv26's own "
                "entry for the full quote)."),
    ),
    'test.tpv31': dict(
        # rule 17 step 5: recorded, not run. requireFaultGeometryResolution
        # admits only dx=500 m today (the gate tier); no 50 m geometry file
        # is shipped in case_input/test.tpv31/ -- a scheduling/storage
        # decision, not a spec gap (planar fault, needs only a finer
        # on_fault_vars grid + material table, not a downloaded surface).
        dx=50.0, term=15.0, nx=4, ny=1, nz=4,
        citation=("TPV31_32_Description_v03, p.9: 'For TPV31, please submit "
                "results using 50 m node spacing on the fault plane.' "
                "'Run the model for times from 0.0 to 15.0 seconds after "
                "nucleation.' ny=1 keeps the fault plane off MPI partitions, "
                "same convention as every other planar-fault full-tier entry "
                "here."),
    ),
    'test.tpv32': dict(
        # Same geometry/term citation as test.tpv31 (identical fault/stress/
        # nucleation/friction, spec Part 2 -- only the 1D velocity structure
        # differs). TPV32's own spec gives a RANGE (25-50 m); 50 m (the
        # coarser end) is recorded here as the scheduling choice (rule 17
        # step 5 is a recording, not a claim this is the spec's preferred
        # value -- see that compset's README).
        dx=50.0, term=15.0, nx=4, ny=1, nz=4,
        citation=("TPV31_32_Description_v03, p.9: 'For TPV32, please select "
                "node spacing on the fault plane in the range of 25 m to "
                "50 m, and submit results for your selected node spacing.' "
                "'Run the model for times from 0.0 to 15.0 seconds after "
                "nucleation.'"),
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
    'test.tpv35': dict(
        # rule 17 step 5: recorded, not run. The case is registered
        # (case_input/test.tpv35, gate dx=500 m at GATE_TERM_S). 100 m is the
        # finest dx the shipped mu_s/tau0 grid admits without interpolation
        # (lib.requireFaultGeometryResolution). ny=1 keeps the fault y-plane
        # off every MPI subdomain boundary, as the gate layout (2,1,2) does;
        # 4x1x2 = 8 ranks for a 60 x 32 x 26 km box at 100 m.
        dx=100., term=18., nx=4, ny=1, nz=2,
        citation=(
            'SCEC TPV35 description (scratch/specs/TPV35_desc.pdf / .txt, '
            'via https://strike.scec.org/cvws/tpv35docs.html -> '
            'download/TPV35_Description_v05.pdf), Part 3 "Running Time, '
            'Node Spacing, and Results": "Run the model for times from 0.0 '
            'to 18.0 seconds after nucleation." / "The recommended '
            'resolution for TPV35 is 100 meters. You may optionally also '
            'submit results for a resolution of 50 meters."'
        ),
    ),
    'test.tpv22': dict(
        # rule 17 step 5: recorded, not run (mission: tpv22/23 campaign,
        # 2026-10-01). ny=1 keeps every fault y-plane off every MPI
        # partition boundary, same reasoning as tpv29/tpv30's full-tier
        # entries -- both our faults sit at distinct, small y-offsets
        # (0 and -1600 m) that an arbitrary y-split could straddle.
        # Geometry note: at this dx the shared x/z mesh belt (see
        # case_input/test.tpv22/tpv22_23_common.py for why it is shared)
        # is 1001 x 401 nodes per fault -- a genuinely large HPC job, same
        # scale class as tpv29's 50 m entry (119 M elements).
        dx=50.0, term=15.0, nx=4, ny=1, nz=4,
        citation=("TPV22_23_Description_v08, p.7 'Running Time, Node "
                "Spacing, and Results': 'Run the model for times from 0.0 "
                "to 15.0 seconds after nucleation.' / 'Please submit "
                "results for two resolutions: Using 100 m node spacing... "
                "Using 50 m node spacing... If you are unable to run the "
                "simulation with 50 m node spacing, then it is OK to "
                "submit just 100 m results.'"),
    ),
    'test.tpv23': dict(
        # Same geometry/duration citation basis as test.tpv22 (shared PDF,
        # shared Part 3/7); only the stepover sign/distance differs
        # (1.0 km compressional vs TPV22's 1.6 km extensional -- Parts 1/2).
        dx=50.0, term=15.0, nx=4, ny=1, nz=4,
        citation=("TPV22_23_Description_v08, p.7 'Running Time, Node "
                "Spacing, and Results' (same page/text as test.tpv22's "
                "entry -- TPV22 and TPV23 share Part 3 verbatim)."),
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
    'test.tpv34': (
        # case_input/test.tpv34 exists (gated at dx=500 m, GATE_TERM_S,
        # 2026-10-06); the shipped CVM-H grid admits only 500 m, so a
        # spec-resolution run needs a fresh extraction with its
        # extract_cvmh_grid.py as well -- another reason not to pick a dx.
        'SCEC TPV34 description (scratch/specs/TPV34_desc.pdf / .txt, via '
        'https://strike.scec.org/cvws/tpv34docs.html -> '
        'download/TPV34_Description_v10.pdf), Part 3 "Running Time, Node '
        'Spacing, and Results": "For TPV34, please select node spacing on '
        'the fault plane in the range of 25 m to 50 m, and submit results '
        'for your selected node spacing." No single recommended value is '
        'given (unlike TPV35\'s "recommended resolution ... is 100 '
        'meters") -- a submitter-chosen range, not a target to rewrite '
        'par.dx to. Excluded rather than invented. Duration is stated '
        '(0 to 20.0 s after nucleation) but withheld here too since dx is '
        'the pairing key for this tier.'
    ),
}
