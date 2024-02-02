# /bin/bash 
# This bash shell script installs necessary packages for 
#   EQdyna to run on a raw Ubuntu 22. 
# It uses MPICH MPI.

apt-get install git vim make mpich
apt-get install libnetcdf-dev libnetcdff-dev 
apt-get install python3 python3-pip
pip install numpy>=1.20 netCDF4 matplotlib xarray