#!/bin/bash
##SBATCH --export=NONE
##SBATCH --get-user-env=L


##
#SBATCH --job-name=TPV2802_50m_a6.slurm
###SBATCH -p mpi
#SBATCH --time=4:30:00
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=48
#SBATCH --mem=360G
#SBATCH --account=132792052222
##SBATCH --partition=short
#SBATCH --output=job.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=d.liu@aggienetwork.com
#

#First Executable Line
#tcsh runcycle.tcsh
date >> log.log
mpirun -np 960 ./run_eqdyna3d
date >> log.log
#
#for ipatch in {39..39}
#do 
#  echo "Patch #:" $ipatch
#  echo "The x"$ipatch" patch started at" >>log.log

#  rm patchnumber.txt
#  rm direction.txt
#  echo $ipatch >>patchnumber.txt
#  for idirection in {1..1}
#	do 
#  	echo $idirection >>direction.txt
#	done
#  mkdir "./x$((ipatch))"
#  mpirun -np 100 ./eqdyna-hyb
#  mv body* frt* fslip* "./x$((ipatch))"
#  
#  echo "y"$ipatch" started at" >>log.log
#  date >> log.log
#  rm direction.txt
#  for idirection in {2..2}
#	do
#  	echo $idirection >>direction.txt
#	done
#  mkdir "./y$((ipatch))"
#  mpirun -np 100 ./eqdyna-hyb
#  mv body* frt* fslip* "./y$((ipatch))"
#  
#  echo "z"$ipatch" patch started at" >>log.log
#  date >> log.log
#  rm direction.txt
#  for idirection in {3..3}
#	do
#  	echo $idirection >>direction.txt
#	done
#  mkdir "./z$((ipatch))"
#  mpirun -np 100 ./eqdyna-hyb
#  mv body* frt* fslip* "./z$((ipatch))"
#done
