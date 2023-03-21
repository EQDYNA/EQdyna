#! user/bin/bash 

module load netcdf
ml

cd src
make
cd ..
mkdir bin
mv src/eqdyna bin

export ESCIROOT=$(pwd)
export PATH=$(pwd)/bin:$PATH
export PATH=$(pwd)/scripts:$PATH

chmod -R 755 scripts
