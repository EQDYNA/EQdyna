#!/bin/bash
# usage: run.sh <case tpv36|tpv37> <dx|gate>
S=<scratch>; T=<worktrees>/cond-screen
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 EQDYNAROOT=$T PATH=$T/bin:$T/scripts:$PATH HWLOC_COMPONENTS="-gl,-opencl,-cuda,-nvml"; unset DISPLAY
d=$S/screen/$1_$2; rm -rf $d; python3 $T/scripts/create.newcase $d test.$1 >/dev/null || exit 3; cd $d
if [ $2 != gate ]; then python3 $S/fine/scripts/set_case.py $d $2 5 1 1 1 || exit 3; else printf '\npar.nx = 1\npar.ny = 1\npar.nz = 1\n' >> user_defined_params.py; fi
python3 case.setup > setup.log 2>&1 || exit 4
timeout 1500 mpirun -np 1 eqdyna --mesh-only > mesh.log 2>&1; echo "rc=$?" ; grep 'mesh-screen' mesh.log
