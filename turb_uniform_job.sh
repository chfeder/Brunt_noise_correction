#!/bin/bash

# Queue type: normal
#PBS -q normal

# 1 job per simulation (21, since 1 resolution * 3 flow types * 7 Mach numbers)
#PBS -J 0-20
#PBS -r y

#PBS -P ek9
#PBS -l walltime=2:00:00
#PBS -l ncpus=1
#PBS -l mem=16GB
#PBS -l storage=scratch/ek9

# Add output files into logs/
#PBS -o logs/

# Working directory
#PBS -l wd

# Job name
#PBS -N turb_uniform

# Join output and error
#PBS -j oe

module load python3/3.12.1
module load openmpi/4.0.7
mkdir -p logs
mpirun -np 1 python3 ./get_turb_data_uniform.py -isim $PBS_ARRAY_INDEX 1>logs/shell.out.$PBS_ARRAY_INDEX 2>&1

# Run this file with the following (since there's a max array size of 10, and min of 2):
# qsub -J 0-9 turb_uniform_job.sh
# qsub -J 10-18 turb_uniform_job.sh
# qsub -J 19-20 turb_uniform_job.sh

# Takes about 50mins to run everything.

# Other commands:
# qstat -u ag9727
# qstat -t "id[]"           = expanded
# nqstat_anu
# ls logs/                  = shell.out.N files
# tail -f logs/shell.out.x  = watch sim x
# qdel "id[]"               = delete a job array