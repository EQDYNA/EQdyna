#!/bin/bash
# usage: launch.sh <tag>   (runs queue.sh on cases_<tag>.txt, detached, cwd = wedge folder)
# ROW149_WEDGE_ROOT: scratch folder for the run directories. ROW149_TREE: a built EQdyna tree carrying the
# row149-force-budget instrumentation patch (bin/eqdyna, scripts/). ROW149_LOCK: lock file shared by heavy runs.
[ -n "$ROW149_WEDGE_ROOT" ] || exit 2
[ -n "$ROW149_TREE" ] || exit 2
W=$ROW149_WEDGE_ROOT
cd "$W" || exit 3
setsid nohup "$W/scripts/queue.sh" "$W/cases_$1.txt" > "$W/queue_$1.log" 2>&1 < /dev/null &
echo $! > "$W/queue_$1.pid"
sleep 2
echo "## launch queue_$1 $(date -Is)" >> "$W/NOTES_wedge.md"
ps -o pid,lstart,args -p "$(cat "$W/queue_$1.pid")" | tee -a "$W/NOTES_wedge.md"
