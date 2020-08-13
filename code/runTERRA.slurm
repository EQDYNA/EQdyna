#!/bin/bash
##SBATCH --export=NONE
##SBATCH --get-user-env=L


##
#SBATCH --job-name=DynamicInversionv5-20200606.slurm
###SBATCH -p mpi
#SBATCH --time=10:00:00
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=20
#SBATCH --mem=50G
#SBATCH --account=122792052368
##SBATCH --partition=short
#SBATCH --output=output.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dunyuliu@tamu.edu
###
echo "srun -l /bin/hostname"
srun -l /bin/hostname

mpirun -np 100 ./eqdyna-hyb
