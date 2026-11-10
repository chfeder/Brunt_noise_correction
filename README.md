# Brunt Noise Correction

This repository contains code for two related projects investigating noise correction in the Brunt Method:

1. **5-Box Project:** Refers to the work described in the paper [*From 2D to 3D: Recovering Turbulent Density Dispersions from Noisy Data*](https://doi.org/10.3847/1538-4357/ae875e).
2. **Uniform Project:** A related investigation using a uniform grid.

## Citation

If you use this code in your research, please cite:

[*From 2D to 3D: Recovering Turbulent Density Dispersions from Noisy Data*](https://ui.adsabs.harvard.edu/abs/2026ApJ..1006..115J/abstract), *The Astrophysical Journal* (2026).

## Workflow Summary for the Uniform Project

1. Simulation data files are soft-linked into `../sims` using `link_sims.py`.
2. Data is processed into a dill file per simulation and placed into `../processed_data` by running `bash/turb_uniform_job.sh`, as explained in the comments of that file. It calls `get_turb_data_uniform.py` for each simulation.
3. Plots are made from the processed data by running `bash/plot_uniform_job.sh`, which executes `plot_turb_data_uniform.py` for each resolution.

## Explanation of Other Files

- `models_driving.txt`: Specifies the simulations to be run.
- `globals_turb_data.py`: Specifies global variables and functions for both projects.
