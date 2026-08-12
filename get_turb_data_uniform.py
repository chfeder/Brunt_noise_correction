#!/usr/bin/env python

from globals_turb_data import *
from link_sims import *
from itertools import product
import flashlib
import subprocess
import dill
import os

# ===================================================================
def get_turb_box(filename, res_resample=1.0):
    # Note res_sample does nothing here currently.
    return flashlib.FlashGG(filename).GetUniformGrid()

# ===================================================================
# ---- generates the model as a string for a given simulation, ie, "N256cM0p2", given the below params
def make_model_name_params(res, flow_type, mach) -> str:
    res = str(res)
    flow_type = str(flow_type)
    mach = float(mach)

    if mach == int(mach):
        mach_str = 'M' + str(int(mach))
    else:
        mach_str = 'M' + str(mach).replace('.', 'p')

    return 'N' + res + flow_type + mach_str

# ===================================================================
def get_output_file_path(project_path, sim_name, noise_exp, res_resample):
    s = f"{project_path}processed_data/{sim_name}"
    if noise_exp != 1.0:
        s += "_noise_exp_"+str(noise_exp)
    if res_resample != 1.0:
        s += "_res_resample_"+str(res_resample)
    s += ".dill"
    return s

# -----------------------------------------------------------

# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    # User changes here:
    resolutions = [256]                         # simulation resolution (default first)
    flow_types = ['s', 'c', 'm']                # Soluniodal, Compressive, Mixture
    mach_nums = [0.1, 0.2, 0.5, 1, 2, 5, 10]    # Mach numbers
    project_path = "/scratch/ek9/ag9727/asc/"   # Path to folder
    min_time = 30                               # Files correponding to times before this will be discarded
    max_time = 1000                             # Likewise, after this will be discarded

    # Test case:
    """
    resolutions = [256]                         # simulation resolution (default first)
    flow_types = ['s']                          # Soluniodal, Compressive, Mixture
    mach_nums = [0.1]                           # Mach numbers
    project_path = "/scratch/ek9/ag9727/asc/"   # Path to folder
    min_time = 30                               # Files correponding to times before this will be discarded
    max_time = 31                               # Likewise, after this will be discarded
    """

    noise_samples = 100 # convergence was checked also with 1000 noise samples, which yields practically the same
    print('Number of noise samples: '+str(noise_samples), color='lightred')

    # Before running the main loop, ensure that all combinations selected would exist
    print("CHECKING EXISTENCE OF MODELS")
    existing_sims = os.listdir(f"{project_path}sims")

    sims = list(product(resolutions, flow_types, mach_nums))
    for res, flow_type, mach in sims:
        sim_name = make_model_name_params(res, flow_type, mach)
        if sim_name in existing_sims:
            print(f"{sim_name} exists.", color='green')
        else:
            print(f"WARNING: {sim_name} DOES NOT EXIST. IT WILL BE SKIPPED.", color='red')


    # loop over all simulations
    timer = cfp.timer("MAIN")

    for res, flow_type, mach in sims:

        # Only try to run the simulation if it exists:
        sim_name = make_model_name_params(res, flow_type, mach)
        if sim_name in existing_sims:
            print(f"===== PROCESSING SIMULATION WITH RES={res}, FLOW={flow_type}, MACH={mach}", color='lightred')
        else:
            print(f"===== SKIPPING SIMULATION WITH RES={res}, FLOW={flow_type}, MACH={mach} (it doesn't exist)", color='red')
            continue

        # Get the names of all timestep files
        timestepped_filenames = [f"{project_path}sims/{sim_name}/Turb_hdf5_plt_cnt_{t:04d}" for t in range(min_time, max_time + 1)]

        # loop over noise exponents
        for noise_exp in [1.0]:

            noise_exp_info = str(noise_exp)
            if noise_exp == 1.0:
                noise_exp_info += ' (white noise)'
            print("===== NOISE EXPONENT = "+noise_exp_info+" =====", color='cyan')

            for res_resample in [1.0]:

                if res_resample != 1.0:
                    print("===== RESOLUTION RESAMPLING FACTOR = "+str(res_resample)+" =====", color='yellow')

                # data collection per box and/or LOS
                data = [[] for _ in range(len(SNR_list))] # The outer axis is SNRs, inner axis is rets

                # Loop over each timestep in the simulation
                for filename in timestepped_filenames:

                    dens = get_turb_box(filename, res_resample=res_resample)
                    print("=== Read dims "+str(dens.shape)+" from file "+filename, color='lightblue')
                    sigma_3D = (dens/dens.mean()).std()

                    # loop over the 3 projection axes (x,y,z)
                    for LOS in LOS_list:
                        print("computing projection along axis "+str(LOS), color='lightblue')

                        proj = dens.sum(axis=LOS) / dens.shape[LOS]

                        # loop over the 7 SNRs
                        for iSNR, SNR in enumerate(SNR_list):

                            ret = noise_iterations(noise_samples, proj, SNR, noise_exp=noise_exp)
                            ret.sigma_3D = np.array([sigma_3D]) # attach sigma_3D

                            data[iSNR].append(ret)

                # compute statistics over noise samples, boxes and LOS'
                stats = []
                for isnr in range(len(SNR_list)):
                    avg, std, p16, p84 = get_stats_object(data[isnr])
                    stats.append(SimpleNamespace(avg=avg, std=std, p16=p16, p84=p84))

                # package and save (dill) all data
                turb_data = {"data": data, "stats": stats}
                dillfile = get_output_file_path(project_path, sim_name, noise_exp, res_resample)
                with open(dillfile, "wb") as f:
                    dill.dump(turb_data, f)
                print("Wrote '"+dillfile+"'.", color='magenta')

    timer.report()