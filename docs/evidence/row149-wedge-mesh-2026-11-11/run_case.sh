#!/bin/bash
# usage: run_case.sh <name> <dx> <nx> <ny> <nz> <np> <timeout_s>
# ROW149_WEDGE_ROOT: scratch folder for the run directories. ROW149_TREE: a built EQdyna tree carrying the
# row149-force-budget instrumentation patch (bin/eqdyna, scripts/). ROW149_LOCK: lock file shared by heavy runs.
[ -n "$ROW149_WEDGE_ROOT" ] || exit 2
[ -n "$ROW149_TREE" ] || exit 2
W=$ROW149_WEDGE_ROOT
T=$ROW149_TREE
name=$1
dx=$2
np=$6
tmo=$7
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export EQDYNAROOT=$T
export PATH=$T/bin:$T/scripts:$PATH
d=$W/$name
if [ -f $d/DONE ]; then echo "$name already DONE, skip"; exit 0; fi
if [ -d $d ]; then mv $d $d.aside.$(date +%s); fi
python3 $T/scripts/create.newcase $d test.tpv12 || exit 3
cd $d || exit 3
python3 $W/scripts/set_case.py $W/params/user_defined_params.py user_defined_params.py $dx $3 $4 $5 || exit 3
./case.setup > case.setup.log 2>&1 || exit 4
echo "start $(date -Is) np=$np" > run.meta
env | grep NUM_THREADS >> run.meta
t0=$(date +%s)
timeout $tmo mpirun -np $np eqdyna > run.log 2>&1 < /dev/null
rc=$?
t1=$(date +%s)
echo "rc=$rc wall=$((t1-t0)) end $(date -Is)" >> run.meta
if [ $rc -eq 0 ]; then touch DONE; fi
cat run.meta
