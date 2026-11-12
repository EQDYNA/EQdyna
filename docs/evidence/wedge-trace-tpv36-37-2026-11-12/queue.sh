#!/bin/bash
# usage: queue.sh <cases_file>; line: name case dx term nx ny nz timeout_s. Needs WT_ROOT WT_TREE WT_LOCK.
[ -n "$WT_LOCK" ] || exit 2
D=$(cd "$(dirname "$0")" && pwd)
while read -r name case dx term nx ny nz tmo; do
    [ -z "$name" ] && continue
    echo "$(date -Is) before $name: $(uptime)"
    load=$(cut -d' ' -f1 /proc/loadavg)
    if ! python3 -c "import sys; sys.exit(0 if float('$load') <= 40 else 1)"; then echo "BLOCKED load $load before $name"; touch "$1.END"; exit 7; fi
    flock "$WT_LOCK" "$D/run_case.sh" $name $case $dx $term $nx $ny $nz $tmo < /dev/null
    echo "$(date -Is) finished $name rc=$?"
done < "$1"
echo "$(date -Is) QUEUE_END $1"; touch "$1.END"
