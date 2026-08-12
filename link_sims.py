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
    if "seed" in model.colnames:
        if model["seed"] == 5887: seed_str = 'S5887'
        if model["seed"] == 1: seed_str = 'S1'
        if model["seed"] == 2: seed_str = 'S2'
        if model["seed"] == 3: seed_str = 'S3'
        if model["seed"] == 4: seed_str = 'S4'
        if model["seed"] == 5: seed_str = 'S5'
        if model["seed"] == 6: seed_str = 'S6'
        if model["seed"] == 7: seed_str = 'S7'
        if model["seed"] == 8: seed_str = 'S8'
        if model["seed"] == 9: seed_str = 'S9'
        if model["seed"] == 140281: seed_str = ''
    return seed_str

def get_mach_str(model):
    mach_str = ""
    if "Mach" in model.colnames:
        if model["Mach"] == 0.1: mach_str = 'M0p1'
        if model["Mach"] == 0.2: mach_str = 'M0p2'
        if model["Mach"] == 0.5: mach_str = 'M0p5'
        if model["Mach"] == 1.0: mach_str = 'M1'
        if model["Mach"] == 2.0: mach_str = 'M2'
        if model["Mach"] == 5.0: mach_str = 'M5'
        if model["Mach"] == 10.0: mach_str = 'M10'
    return mach_str

def get_aviro_str(model):
    aviro_str = ""
    if "aviro" in model.colnames:
        if model["aviro"] == 0.5: aviro_str = 'A0p5'
        if model["aviro"] == 1.0: aviro_str = 'A1'
        if model["aviro"] == 2.5: aviro_str = 'A2p5'
        if model["aviro"] == 5.0: aviro_str = 'A5'
        if model["aviro"] == 10.0: aviro_str = 'A10'
        if model["aviro"] == 20.0: aviro_str = 'A20'
    else:
        aviro_str = ""
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