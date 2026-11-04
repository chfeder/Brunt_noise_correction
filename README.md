# Brunt Noise Correction

## Workflow Summary for the Uniform Project

1. Soft-link simulations into ../sims using `link_sims.py`
2. Process data into a dill file per simulation into ../processed_data by running `bash/turb_uniform_job.sh` as explained in the comments for that file. It calls `get_turb_data_uniform.py` for each simulation.
3. Convert the processed data into simulations by running `bash/plot_uniform_job.sh`, which executes the code in `plot_turb_data_uniform.py` for each resolution.

## Explaination of Other Files
    - `models_driving.txt`:     Specifies the simulations to be run.
    - `globals_turb_data.py`:   Specifies global variables and functions for both projects.


# Variable Renaming


In `plot_turb_data_uniform.py`, I've renamed many of the modelling variables in the labels but haven't
updated the code to reflect that changes, which may be confusing to read. Here's a table showing old
vs new paramater names

| Old | New | Role |
|---|---|---|
| $\mathrm{SNR_{crit}}$ | $\mathrm{SNR_{crit}}$ | unchanged |
| $\eta$ | $a$ | input to $\mathrm{SNR_{crit}}$ |
| $\lambda_{\text{old}}$ | $b$ | input to $\mathrm{SNR_{crit}}$ |
| $PL_EXP = \alpha$ | $\lambda$ | steepness (the exponent) |
| $\alpha_0$ | $c$ | input to $\lambda$ |
| $\alpha_1$ | $d$ | input to $\lambda$ |