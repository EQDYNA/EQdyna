#!/bin/bash
# usage: queue.sh <cases_file>; each line: name dx nx ny nz np timeout_s
# ROW149_WEDGE_ROOT: scratch folder for the run directories. ROW149_TREE: a built EQdyna tree carrying the
# row149-force-budget instrumentation patch (bin/eqdyna, scripts/). ROW149_LOCK: lock file shared by heavy runs.
[ -n "$ROW149_WEDGE_ROOT" ] || exit 2
[ -n "$ROW149_TREE" ] || exit 2
W=$ROW149_WEDGE_ROOT
[ -n "$ROW149_LOCK" ] || exit 2
L=$ROW149_LOCK
while read -r name dx nx ny nz np tmo; do
    [ -z "$name" ] && continue
    echo "$(date -Is) before $name: $(uptime) nproc $(nproc)"
    load=$(cut -d' ' -f1 /proc/loadavg)
    if ! python3 -c "import sys; sys.exit(0 if float('$load') <= 48 else 1)"; then echo "BLOCKED load $load before $name" | tee "$W/QUEUE_BLOCKED"; touch "$1.END"; exit 7; fi
    flock "$L" "$W/scripts/run_case.sh" "$name" "$dx" "$nx" "$ny" "$nz" "$np" "$tmo" < /dev/null
    echo "$(date -Is) finished $name rc=$?"
done < "$1"
echo "$(date -Is) QUEUE_END $1"
touch "$1.END"
