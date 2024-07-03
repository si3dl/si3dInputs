#!/bin/bash -l
# NOTE the -l flag!

# If you need any help, please email help@cse.ucdavis.edu
 
# Name of the job 
#SBATCH -J llq1
# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o llq1-%j-%j.output
#SBATCH -e llq1-%j.output

# no -n here, the user is expected to provide that on the command line.

# The useful part of your job goes below

# run one thread for each one the user asks the queue for
# hostname is just for debugging
hostname
export OMP_NUM_THREADS=$SLURM_NTASKS
module load benchmarks intel

# The main job executable to run: note the use of srun before it
time srun psi3d
