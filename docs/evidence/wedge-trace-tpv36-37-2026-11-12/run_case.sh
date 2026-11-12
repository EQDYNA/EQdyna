#!/bin/bash
# usage: run_case.sh <name> <tpv36|tpv37> <dx> <term_s> <nx> <ny> <nz> <timeout_s>
# WT_ROOT: scratch folder for run dirs. WT_TREE: built EQdyna tree (bin/eqdyna, scripts/). Ranks = nx*ny*nz.
[ -n "$WT_ROOT" ] && [ -n "$WT_TREE" ] || exit 2
W=$WT_ROOT; T=$WT_TREE
name=$1; case=$2; dx=$3; term=$4; nx=$5; ny=$6; nz=$7; tmo=$8
np=$((nx*ny*nz))
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export EQDYNAROOT=$T PATH=$T/bin:$T/scripts:$PATH
d=$W/$name
if [ -f $d/DONE ]; then echo "$name already DONE, skip"; exit 0; fi
if [ -d $d ]; then mv $d $d.aside.$(date +%s); fi
python3 $T/scripts/create.newcase $d test.$case || exit 3
cd $d || exit 3
python3 "$(dirname "$0")/set_case.py" $d $dx $term $nx $ny $nz || exit 3
./case.setup > case.setup.log 2>&1 || exit 4
echo "start $(date -Is) np=$np" > run.meta
env | grep NUM_THREADS >> run.meta
t0=$(date +%s)
timeout $tmo mpirun -np $np eqdyna > run.log 2>&1 < /dev/null
rc=$?
echo "rc=$rc wall=$(( $(date +%s)-t0 )) end $(date -Is)" >> run.meta
[ $rc -eq 0 ] && touch DONE
cat run.meta
