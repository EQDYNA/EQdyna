#!/bin/bash
#SBATCH -J TPV2802_50_b            # job name
#SBATCH -o mpiout.o%j        # output and error file name (%j expands to jobID)
#SBATCH -N 20                # number of nodes requested
#SBATCH -n 960              # total number of mpi tasks requested
#SBATCH -p skx-normal      # queue (partition) --small (2), normal(3-512), development, etc.
#SBATCH -t 5:00:00         # run time (hh:mm:ss) - 1.5 hours
#SBATCH -A High-frenquency-eart
#SBATCH --mail-user=dliu@ig.utexas.edu
#SBATCH --mail-type=begin   # email me when the job starts
#SBATCH --mail-type=end     # email me when the job finishes

date >> log.log
ibrun run_eqdyna3d
date >> log.log
