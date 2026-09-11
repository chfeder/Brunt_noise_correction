#!/bin/bash

# Queue type: normal
#PBS -q normal

# 1 job per resolution (128, 256)
#PBS -r y

#PBS -P ek9
#PBS -l walltime=00:20:00
#PBS -l ncpus=1
#PBS -l mem=128GB
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
mpirun -np 1 python3 ./plot_turb_data_uniform.py -ires 0 1>logs/plot.out.0 2>&1

# Run this file with:
# qsub ./plot_uniform_job.sh

# usually takes about 4 mins