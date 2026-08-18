#!/bin/bash

# Queue type: normal
#PBS -q normal


#PBS -P ek9
#PBS -l walltime=00:10:00
#PBS -l ncpus=1
#PBS -l mem=64GB
#PBS -l storage=scratch/ek9

# Add output files into logs/
#PBS -o logs/

# Working directory
#PBS -l wd

# Job name
#PBS -N plot_turb

# Join output and error
#PBS -j oe

module load python3/3.12.1
module load openmpi/4.0.7
mkdir -p logs
mkdir -p ../figs
mpirun -np 1 python3 ./plot_turb_data_uniform.py 1>logs/plot.out 2>&1

# Run this file with the following
# qsub plot_uniform_job.sh