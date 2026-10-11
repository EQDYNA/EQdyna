#!/bin/bash
# usage: runall.sh CASE TERM VARIANTS...   (variants: base off on)
set -e
WT=${WT:?set WT to the EQdyna checkout}
S=${S:?set S to a scratch dir holding bin/eqdyna.base and runs/}
CASE=$1; TERM_S=$2; NAME=${NAME:-$CASE}; shift 2
export EQDYNAROOT=$WT PATH=$WT/bin:$WT/scripts:$PATH OMP_NUM_THREADS=1 PYTHONPATH=$WT/src/python
cd $S/runs
rm -rf $NAME.setup
nice create.newcase $NAME.setup $CASE > $NAME.setup.log 2>&1
printf '\npar.term = %s\npar.outputFaultTimeSeries = 1\n' $TERM_S >> $NAME.setup/user_defined_params.py; [ -n "$EXTRA" ] && printf "\n%s\n" "$EXTRA" >> $NAME.setup/user_defined_params.py
(cd $NAME.setup && nice ./case.setup >> ../$NAME.setup.log 2>&1)
for v in "$@"; do
  d=$NAME.$v; rm -rf $d; cp -r $NAME.setup $d
  if [ $v != on ]; then sed -i '$ s/^1$/0/' $d/bGlobal.txt; fi
  bin=$WT/bin/eqdyna; [ $v = base ] && bin=$S/bin/eqdyna.base
  (cd $d && /usr/bin/time -f "%e s" nice mpirun -np 4 $bin > run.log 2>&1 && nice python3 plotRuptureDynamics > plot.log 2>&1) || { echo "$d FAILED"; tail -20 $d/run.log; exit 1; }
  echo "$d done $(tail -1 $d/run.log)"
done
