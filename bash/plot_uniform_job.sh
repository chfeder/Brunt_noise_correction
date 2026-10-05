#!/bin/bash

# Normal queue
#PBS -q normal

# 2 jobs - 0 for "SMC Data", and 1 for "Zeta Data" (explained in pycode)
#PBS -J 0-1

# Rerunnable
#PBS -r y

#PBS -P ek9
#PBS -l walltime=00:15:00
#PBS -l ncpus=1
#PBS -l mem=192GB
#PBS -l storage=scratch/ek9

# Add output files into logs/
#PBS -o ../logs/

# Working directory
#PBS -l wd

# Job name
#PBS -N plot_turb

# Join output and error
#PBS -j oe

module load python3/3.12.1
module load openmpi/4.0.7
mkdir -p ../logs
mkdir -p ../figs
mkdir -p ../figs/main
mkdir -p ../figs/extra

# Each sub-job gets its own index (0 or 1) from PBS
IRES=${PBS_ARRAY_INDEX}
mpirun -np 1 python3 ../plot_turb_data_uniform.py -ires ${IRES} 1>../logs/plot.out.${IRES} 2>&1

# Run this file with:
# qsub ./plot_uniform_job.sh

# Schedule with:
# qsub -a $(date -d "+3 hours" +%Y%m%d%H%M) plot_uniform_job.sh

# Watch execution:
# tail -f ../logs/plot.out.0
# tail -f ../logs/plot.out.1

# Check sub-jobs:
# qstat -t

# usually takes about 4 mins