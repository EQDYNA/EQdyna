#! /bin/bash 

module load netcdf
ml

echo $TACC_NETCDF_INC
echo $TACC_NETCDF_LIB
export FC=mpiifort

cd src
make
cd ..
mkdir bin
mv src/eqdyna bin

export ESCIROOT=$(pwd)
export PATH=$(pwd)/bin:$PATH
export PATH=$(pwd)/scripts:$PATH

chmod -R 755 scripts
