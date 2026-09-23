#! /bin/sh
# Regenerate the SCEC cross-version comparison proof set with scec_compare.py.
#
# Every figure here comes from ONE reusable script; nothing is hand-edited and
# no figure exists without this driver being able to remake it.
#
#   sh scripts/figures/make_scec_comparison_set.sh [OUTDIR]
#
# EQDYNA_MAIN must point at a checkout that HAS scec_archive/ and the run
# directories -- scec_archive/ is gitignored (.gitignore:39), so it lives only
# in the primary checkout and is never copied into a worktree. It is opened
# READ-ONLY by every command below.
set -e
MAIN=${EQDYNA_MAIN:-/home/utig5/dliu/EQdyna}
A=$MAIN/scec_archive
OUT=${1:-scratch/scec_compare_proof}
S=$(dirname "$0")/scec_compare.py
mkdir -p "$OUT"

# --- TPV30: our 500 m run against ALL THREE archived resolutions at once ----
# The convergence picture: 100 m, 50 m and 25 m submissions of the same
# benchmark overlaid on one axes with our run.
python3 "$S" --plot cplot --models \
    "$MAIN/test/test.tpv30" \
    "$A/tpv30/eqdyna-v3.1-100m-2015" \
    "$A/tpv30/eqdyna-v3.1-50m-2015" \
    "$A/tpv30/eqdyna-v3.1-25m-2015" \
    --out "$OUT/tpv30_cplot.png"
python3 "$S" --plot ts-fault --models \
    "$MAIN/test/test.tpv30" \
    "$A/tpv30/eqdyna-v3.1-100m-2015" \
    "$A/tpv30/eqdyna-v3.1-50m-2015" \
    --out "$OUT/tpv30_ts_fault.png"
python3 "$S" --plot ts-body --models \
    "$MAIN/test/test.tpv30" \
    "$A/tpv30/eqdyna-v3.1-100m-2015" \
    "$A/tpv30/eqdyna-v3.1-50m-2015" \
    --out "$OUT/tpv30_ts_body.png"

# --- TPV104: different fault extent, 9-station set, rate-and-state ----------
python3 "$S" --plot cplot --models \
    "$MAIN/test.reference.results/test.tpv104" \
    "$A/tpv104/eqdyna3d-v4.1-100m-2016" \
    --out "$OUT/tpv104_cplot.png"

# --- TPV105-3d: different extent again, 12-13 station set ------------------
python3 "$S" --plot cplot --models \
    "$MAIN/scratch/tpv1053d.diffcheck/test.tpv1053d" \
    "$A/tpv105-3d/eqdyna3d-v5.1.0-100m-2020" \
    "$A/tpv105-3d/eqdyna3d-v5.1.0-200m-2020" \
    --out "$OUT/tpv1053d_cplot.png"
python3 "$S" --plot ts-fault --models \
    "$MAIN/scratch/tpv1053d.diffcheck/test.tpv1053d" \
    "$A/tpv105-3d/eqdyna3d-v5.1.0-100m-2020" \
    "$A/tpv105-3d/eqdyna3d-v5.1.0-200m-2020" \
    --out "$OUT/tpv1053d_ts_fault.png"
python3 "$S" --plot ts-body --models \
    "$MAIN/scratch/tpv1053d.diffcheck/test.tpv1053d" \
    "$A/tpv105-3d/eqdyna3d-v5.1.0-100m-2020" \
    --out "$OUT/tpv1053d_ts_body.png"

# --- TPV36: DIPPING fault, and the 2024-era 99999 s sentinel on both sides --
# Exercises the down-dip derivation (hypot(y,z), decided from the node cloud)
# and the per-file sentinel detection against a different fill value.
python3 "$S" --plot cplot --models \
    "$MAIN/test.reference.results/test.tpv36" \
    "$A/tpv36/eqdyna-v5.3.3-50m-2024" \
    --out "$OUT/tpv36_cplot.png"

echo "regenerated the proof set in $OUT"
