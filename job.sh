#!/bin/bash
#PBS -P ek9
#PBS -q normal
#PBS -l walltime=6:00:00
#PBS -l ncpus=8
#PBS -l mem=3GB
#PBS -l storage=scratch/ek9+gdata/ek9
#PBS -l wd
#PBS -N snr_iterations
#PBS -j oe
#PBS -m bea
#PBS -M ljimenezvela@fsu.edu

source /g/data/ek9/yt-conda/bin/activate
python main_snr_getdata.py 
