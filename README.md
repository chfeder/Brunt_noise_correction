# Brunt Noise Correction

## Workflow Summary for the Uniform Project

1. Soft-link simulations into ../sims using `link_sims.py`
2. Process data into a dill file per simulation into ../processed_data by running `turb_uniform_job.sh` as explained in the comments for that file. It calls `get_turb_data_uniform.py` for each simulation.
3. Convert the processed data into simulations by running `plot_uniform_job.sh`, which executes the code in `plot_turb_data_uniform.py` for each resolution.

## Explaination of Other Files
    - `models_driving.txt`:     Specifies the simulations to be run.
    - `globals_turb_data.py`:   Specifies global variables and functions for both projects.