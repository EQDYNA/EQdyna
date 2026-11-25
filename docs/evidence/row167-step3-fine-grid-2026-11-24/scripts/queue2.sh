#!/bin/bash
export HWLOC_COMPONENTS="-gl,-opencl,-cuda,-nvml"; unset DISPLAY
# serial queue, waits (never times out) for load<=40 and for the PR-head e2e + 500 m run to end before each case
S=<scratch>; Q=<scratch>/fine; EW=<worktrees>
waitfree() { while kill -0 $(cat $S/e2e_pr284.pid) 2>/dev/null || kill -0 $(cat $S/tpv12wedge/q500.pid) 2>/dev/null || ! python3 $Q/idle_ok.py; do sleep 60; done; }
waitfree; echo "$(date -Is) start tpv12 wedge 200 m: $(uptime)"
ROW149_WEDGE_ROOT=$S/tpv12wedge ROW149_TREE=$EW/cond-e flock $Q/lock $S/tpv12wedge/scripts/run_case.sh w12_dx200_master 200 2 2 1 4 5400 < /dev/null
for spec in "250 tpv36 13" "250 tpv37 9" "200 tpv36 13" "200 tpv37 9"; do set -- $spec; dx=$1; c=$2; term=$3
  for side in post pre; do
    tree=$EW/cond-pr; [ $side = pre ] && tree=$EW/cond-pre
    waitfree; echo "$(date -Is) start ${side}_${c}_dx$dx: $(uptime)"
    WT_ROOT=$Q/runs WT_TREE=$tree flock $Q/lock $Q/scripts/run_case.sh ${side}_${c}_dx$dx $c $dx $term 4 1 4 40000 < /dev/null
    echo "$(date -Is) finished ${side}_${c}_dx$dx rc=$?"
  done
done
echo "$(date -Is) QUEUE_END"
