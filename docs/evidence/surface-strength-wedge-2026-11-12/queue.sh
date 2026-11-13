#!/bin/bash
# usage: queue.sh <cases_file>; line: name kind dx nx ny nz term timeout_s. Needs SS_ROOT SS_TREE SS_WPARAMS SS_LOCK.
# Before each run: wait (<= 10 min) until 1-min load + this run's ranks <= 60 on the 64-core box, else stop BLOCKED.
[ -n "$SS_LOCK" ] || exit 2
H=$(dirname "$(readlink -f "$0")")
while read -r name kind dx nx ny nz term tmo; do
    [ -z "$name" ] && continue
    np=$((nx*ny*nz)); ok=0
    for i in $(seq 1 20); do
        load=$(cut -d' ' -f1 /proc/loadavg)
        if python3 -c "import sys; sys.exit(0 if float('$load')+$np <= 60 else 1)"; then ok=1; break; fi
        sleep 30
    done
    echo "$(date -Is) before $name: load $load np $np nproc $(nproc)"
    [ $ok -eq 1 ] || { echo "BLOCKED load $load before $name"; touch "$1.END"; exit 7; }
    flock "$SS_LOCK" "$H/run_case.sh" $name $kind $dx $nx $ny $nz $term $tmo < /dev/null
    echo "$(date -Is) finished $name rc=$?"
done < "$1"
echo "$(date -Is) QUEUE_END $1"; touch "$1.END"
