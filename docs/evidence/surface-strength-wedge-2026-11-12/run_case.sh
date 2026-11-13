#!/bin/bash
# usage: run_case.sh <name> <kind> <dx> <nx> <ny> <nz> <term> <timeout_s>
# SS_ROOT scratch run folder, SS_TREE built tree (bin/eqdyna, scripts/), SS_WPARAMS row149 wedge params file.
[ -n "$SS_ROOT" ] && [ -n "$SS_TREE" ] && [ -n "$SS_WPARAMS" ] || exit 2
name=$1; kind=$2; dx=$3; nx=$4; ny=$5; nz=$6; term=$7; tmo=$8; np=$((nx*ny*nz))
H=$(dirname "$(readlink -f "$0")")
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export EQDYNAROOT=$SS_TREE PATH=$SS_TREE/bin:$SS_TREE/scripts:$PATH
case $kind in w12*|i12*) comp=test.tpv12;; t36*) comp=test.tpv36;; *) exit 2;; esac
d=$SS_ROOT/$name
if [ -f $d/DONE ]; then echo "$name already DONE, skip"; exit 0; fi
if [ -d $d ]; then mv $d $d.aside.$(date +%s); fi
python3 $SS_TREE/scripts/create.newcase $d $comp > /dev/null || exit 3
cd $d || exit 3
python3 "$H/set_case.py" $kind $d $dx $nx $ny $nz $term $SS_WPARAMS > set_case.log 2>&1 || { cat set_case.log; exit 3; }
./case.setup > case.setup.log 2>&1 || exit 4
{ echo "start $(date -Is) np=$np kind=$kind"; env | grep NUM_THREADS; cat set_case.log; } > run.meta
t0=$(date +%s)
timeout $tmo mpirun -np $np eqdyna > run.log 2>&1 < /dev/null &
pid=$!; echo $pid > run.pid
# watchdog: a warning line repeated > 50 times stops the run (by this job's own PID, args confirmed first)
while kill -0 $pid 2>/dev/null; do
    sleep 20
    rep=$(grep -i 'warn' run.log | sort | uniq -c | sort -rn | head -1 | awk '{print $1+0}')
    if [ "${rep:-0}" -gt 50 ] && ps -o pid,args -p $pid | grep -q "timeout $tmo mpirun"; then
        echo "KILL repeated warning x$rep" >> run.meta; kill -TERM $pid; fi
done
wait $pid; rc=$?
echo "rc=$rc wall=$(( $(date +%s)-t0 )) end $(date -Is)" >> run.meta
[ $rc -eq 0 ] && touch DONE
echo "$name,$kind,$dx,$term,$np,$(( $(date +%s)-t0 )),$rc,$(date -Is)" >> $SS_ROOT/runs.csv
tail -1 run.meta
