#! /bin/bash 

# This script will perform tests on default test cases.
# Currently, it includes tpv104 and tpv1053d with smaller domains and coarse element sizes.
#

MPIRUN=mpirun.mpich

rm -rf test
mkdir test
cd test

# create and test case tpv104 with 4 cores
create.newcase tpv104 test-tpv104
cd tpv104 
./case.setup 
$MPIRUN -np 4 eqdyna
python3 plotRuptureDynamics
cd ..
# need to add a script to calculate errors

create.newcase tpv1053d test-tpv1053d
cd tpv1053d
./case.setup 
$MPIRUN -np 4 eqdyna
python3 plotRuptureDynamics
cd ..

create.newcase meng2023a test-meng2023a
cd meng2023a
./case.setup 
$MPIRUN -np 4 eqdyna
python3 plotRuptureDynamics
cd ..

create.newcase meng2023cb test-meng2023cb
cd meng2023cb
./case.setup 
$MPIRUN -np 4 eqdyna
python3 plotRuptureDynamics
cd ..

cd ..



