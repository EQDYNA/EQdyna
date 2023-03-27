#! /bin/bash 

module load netcdf
ml

export TACC_NETCDF_INC=/usr/include
export TACC_NETCDF_LIB=/usr/lib/x86_64-linux-gnu
export FC=mpif90

cd src
make
cd ..
mkdir bin
mv src/eqdyna bin

export EQDYNAROOT=$(pwd)
export PATH=$(pwd)/bin:$PATH
export PATH=$(pwd)/scripts:$PATH

chmod -R 755 scripts
