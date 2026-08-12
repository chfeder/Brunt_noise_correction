#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Christoph Federrath, 2021-2026

from cfpack import stop, print
import cfpack as cfp
import cfpack.constants as myconst
import numpy as np
import os, shutil
import glob
import subprocess
import argparse
import ics
import flashlib as fl


# get the initial conditions
def get_ics(model):
    if make_model_name(model)[0] == "r": # star formation (resim) model
        L = 6.172e+18 # size of box
        cs = 0.2e5 # sound speed
        sigma_v = cs * model["Mach"] # turbulent velocity dispersion
        rho = 5.0 * sigma_v**2 / (6.0 * myconst.g_n * L**2 * model["aviro"]) # get density from Mach and virial parameter
        ics_derived = ics.get_ics(resolution=model["resolution"], L=L, domain="box", rho=rho, sigma_v=sigma_v, quiet=True)
    else: # driving model
        L = 1.0 # size of box
        cs = 1.0 # sound speed
        sigma_v = cs * model["Mach"] # turbulent velocity dispersion
        rho = 1.0 # density
        ics_derived = None
    ics_ret = {"L": L}
    ics_ret["cs"] = cs
    ics_ret["sigma_v"] = sigma_v
    ics_ret["t_turb"] = L / (2.0*sigma_v)
    ics_ret["rho"] = rho
    if ics_derived is not None:
        for k in ics_derived.keys(): # attach derived ics properties
            ics_ret[k] = ics_derived[k]
    return ics_ret

def get_driving_str(model):
    driving_str = ""
    if "driving" in model.colnames:
        if model["driving"] == "sol" : driving_str = 's'
        if model["driving"] == "mix" : driving_str = 'm'
        if model["driving"] == "comp": driving_str = 'c'
    return driving_str

def get_seed_str(model):
    seed_str = ""
    if ("seed" in model.colnames) and (model["seed"] != 140281):
        seed_str = f"S{model["seed"]}"
    return seed_str

def get_mach_str(model):
    mach_str = ""
    if "Mach" in model.colnames:
        mach = model["Mach"]
        if mach == int(mach):
            mach_str = 'M' + str(int(mach))
        else:
            mach_str = 'M' + str(mach).replace('.', 'p')
    return mach_str

def get_aviro_str(model):
    aviro_str = ""
    if "aviro" in model.colnames:
        aviro = model["aviro"]
        if aviro == int(aviro):
            aviro_str = 'A' + str(int(aviro))
        else:
            aviro_str = 'A' + str(aviro).replace('.', 'p')
    return aviro_str

def make_model_name(model):
    if get_aviro_str(model) == "":
        pre_str = "" # driving model
    else:
        pre_str = "r" # star formation (resim) model
    model_name = pre_str+'N'+str(model["resolution"])+get_driving_str(model)+get_mach_str(model)+\
                 get_aviro_str(model)+get_seed_str(model)
    return model_name


# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    # parse command line arguments
    parser = argparse.ArgumentParser(description='Prepare FLASH simulation(s).')
    parser.add_argument("-model_file", type=str, default="models_driving.txt", help="Model filename. Default '%(default)s'")
    parser.add_argument("-submit_job", action='store_true', default=False, help="After preparing files, submit job. Default '%(default)s'")
    parser.add_argument("-overwrite", action='store_true', default=False, help="Overwrite model(s). Default '%(default)s'")
    parser.add_argument("-num_t_turb", type=float, default=10.0, help="Number of turbulent turnover times. Default '%(default)s'")
    parser.add_argument("-adjust_amplitude", action='store_true', default=False, help="Adjust turbulence driving amplitude and rerun. Default '%(default)s'")
    args = parser.parse_args()

    # read simulation models file
    models = cfp.read_ascii(args.model_file)

    # loop over all models
    for i in range(len(models)):

        m = models[i]
        model_name = make_model_name(m)

        print("=== Working on model "+model_name+" ===", color='green')

        source_path = "/scratch/ek9/ccf100/ari/"
        destination_path = "/scratch/ek9/ag9727/asc/sims/"

        cfp.run_shell_command("mkdir " + destination_path + model_name)
        cfp.run_shell_command("ln -sf " + source_path + model_name + "/* " + destination_path + model_name)