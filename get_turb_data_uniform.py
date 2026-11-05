#!/usr/bin/env python

"""
This file contains the code which creates the processed data files for the uniform project given the
simulations as input.

The data files are stored as .dill files, with the following structure:

{"data": data, "stats": stats}

Each data[iSNR][j] is a ret object, where each iSNR is a SNR, and each j is a projection.
The ret object holds the data for all noise_samples iterations. Each ret object contains
the data for that projection.

Each stats[iSNR] is a list, [avg, std, p16, p64], where each entry is a ret object,
containing the data aggregated over the associated statistical metric.
"""

from globals_turb_data import *
from link_sims import *
from itertools import product
import flashlib
import subprocess
import dill
import os
import argparse

# ===================================================================
def get_turb_box(filename, res_resample=1.0):
    # Note res_sample does nothing here currently.
    return flashlib.FlashGG(filename).GetUniformGrid()

# ===================================================================
def get_output_file_path(project_path, sim_name, noise_exp, res_resample):
    s = f"{project_path}processed_data/{sim_name}"
    if noise_exp != 1.0:
        s += "_noise_exp_"+str(noise_exp)
    if res_resample != 1.0:
        s += "_res_resample_"+str(res_resample)
    s += ".dill"
    return s

# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    # Parse the -isim argument (if it exists).
    # This is used for Parallelising, by assigning each sim to a different core.
    parser = argparse.ArgumentParser()
    parser.add_argument("-isim", type=int, default=None, help="index of single sim to process (for job arrays)")
    args = parser.parse_args()

    noise_samples = 100 # convergence was checked also with 1000 noise samples, which yields practically the same
    print('Number of noise samples: '+str(noise_samples), color='lightred')

    # Before running the main loop, ensure that all combinations selected would exist
    print("CHECKING EXISTENCE OF MODELS")
    existing_sims = os.listdir(f"{project_path}sims")

    sims = list(product(resolutions, flow_types, mach_nums))
    # Reduce the list of simulations to do to just the 1 per core if isim is used
    if args.isim is not None:
        sims = [sims[args.isim]]

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
                data = [[] for _ in range(len(SNR_list_uniform))] # The outer axis is SNRs, inner axis is rets

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
                        for iSNR, SNR in enumerate(SNR_list_uniform):

                            # Add the noise randomly for a total of noise_samples times, then store it
                            ret = noise_iterations(noise_samples, proj, SNR, noise_exp=noise_exp)
                            ret.sigma_3D = np.array([sigma_3D]) # attach sigma_3D
                            data[iSNR].append(ret)

                # compute statistics
                stats = []
                for iSNR in range(len(SNR_list_uniform)):
                    avg, std, p16, p84 = get_stats_object(data[iSNR]) # data[iSNR] is all the rets per projection for iSNR
                    stats.append(SimpleNamespace(avg=avg, std=std, p16=p16, p84=p84))

                # package and save (dill) all data
                turb_data = {"data": data, "stats": stats}
                dillfile = get_output_file_path(project_path, sim_name, noise_exp, res_resample)
                with open(dillfile, "wb") as f:
                    dill.dump(turb_data, f)
                print("Wrote '"+dillfile+"'.", color='magenta')


    timer.report()
    print("End.")