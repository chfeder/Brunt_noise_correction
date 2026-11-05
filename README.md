# Brunt Noise Correction

This repository contains code for two related projects surrounding noise correction in the Brunt Method:
1. The '5-Box Project' (refering to the work described in the paper "From 2D to 3D: Recovering Turbulent Density Dispersions from Noisy Data", https://doi.org/10.3847/1538-4357/ae875e)
2. The 'Uniform Project' (a related investigation, using a uniform grid)

## Workflow Summary for the Uniform Project

1. Simulation data files are soft-linked into ../sims using `link_sims.py`
2. Data is processed into a dill file per simulation into `../processed_data by` running `bash/turb_uniform_job.sh` as explained in the comments for that file. It calls `get_turb_data_uniform.py` for each simulation.
3. Plots are made from the processed data by running `bash/plot_uniform_job.sh`, which executes the code in `plot_turb_data_uniform.py` for each resolution.

## Explaination of Other Files
    - `models_driving.txt`:     Specifies the simulations to be run.
    - `globals_turb_data.py`:   Specifies global variables and functions for both projects.