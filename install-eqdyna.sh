#! user/bin/bash 

module load netcdf
ml

cd src
make
cd ..
mkdir bin
mv src/eqdyna bin

export PATH=$(pwd)/bin:$PATH
export PATH=$(pwd)/scripts:$PATH
export ECCIROOT=$(pwd)
