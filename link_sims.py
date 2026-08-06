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

        # New code:

        source_path = "/scratch/ek9/ccf100/ari/"
        destination_path = "/scratch/ek9/ag9727/asc/sims/"

        cfp.run_shell_command("mkdir " + destination_path + model_name)
        cfp.run_shell_command("ln -sf " + source_path + model_name + "/* " + destination_path + model_name)

        # Old code:
        """
        star_formation_model = False
        if model_name[0] == "r": # star formation (resim) model
            star_formation_model = True

        path = '../'+model_name
        if not os.path.isdir(path):
            print('creating dir '+path)
            os.mkdir(path)
        else:
            if not args.adjust_amplitude:
                # check if this model has run already
                shellouts = sorted(glob.glob(path+"/shell.out*"))
                if len(shellouts) > 0:
                    print("!!! Model '"+path+"' has shell outputs --- Skipping !!!", color='lightred')
                    continue
                if args.overwrite:
                    print("!!! Deleting model '"+path+"' !!!", color='lightred')
                    cmd = "rm -rf "+path+"/* "+path+"/.*"
                    print("now running '"+cmd+"'", color='lightred')
                    sp_result = subprocess.run(cmd, shell=True)
                else:
                    key_input = input("Model exisits! This will overwrite lines in flash.par, job.sh, ics.sh, "+
                                    "and forcing_generator.inp. Press 'p' to proceed or any other key to skip...")
                    if key_input != 'p':
                        print("Skipping this model.")
                        continue

        # copy driving template files
        if not args.adjust_amplitude:
            if star_formation_model:
                template_folder = 'template_starform'
            else:
                template_folder = 'template_driving'
            print("copying template files from '../"+template_folder+"'...")
            for f in glob.glob('../'+template_folder+'/*'):
                shutil.copy(f, path)

            # determine initial conditions
            ic = get_ics(m)
            if star_formation_model:
                print("Virial parameter (recomputed) = "+str(cfp.round(ic["virial_parameter"],3))) # confirm virial parameter

        # edit turbulence generator input file
        driv_gen_file = path+'/TurbGen.par'
        if not args.adjust_amplitude:
            new_line = "L                  = "+cfp.eform(ic["L"],6,print_leading_plus=True)
            cfp.replace_line_in_file(driv_gen_file, 'L', new_line)
            new_line = "velocity           = "+str(ic["sigma_v"])+"        ! Target velocity dispersion"
            cfp.replace_line_in_file(driv_gen_file, 'velocity', new_line)
            # st_solweight (sol, mix, comp)
            if m["driving"] == "sol" : solw = 1.0
            if m["driving"] == "mix" : solw = 0.5
            if m["driving"] == "comp": solw = 0.0
            new_line = "sol_weight         = "+str(solw)+"       ! 1.0: solenoidal driving, 0.0: compressive driving, 0.5: natural mixture"
            cfp.replace_line_in_file(driv_gen_file, 'sol_weight', new_line)
            new_line = "random_seed        = "+str(m["seed"])+"    ! Random seed for forcing sequence"
            cfp.replace_line_in_file(driv_gen_file, 'random_seed', new_line)
            new_line = "ampl_factor        = "+str(1.0)+"    ! Amplitude factor"
            cfp.replace_line_in_file(driv_gen_file, 'ampl_factor', new_line)
            new_line = "ampl_auto_adjust   = "+str(1)+"    ! Amplitude auto adjustment"
            cfp.replace_line_in_file(driv_gen_file, 'ampl_auto_adjust', new_line)
        else:
            datfile = path+"/TurbGen.dat"
            if not os.path.isfile(datfile): print("Need to run auto adjustment of amplitude first", error=True)
            TGdatfile = fl.datfile(datfile)
            ampl_x_mean = TGdatfile.get_moments("#03_ampl_factor_x", xcol=1, xs=3.0, xe=args.num_t_turb)['mean']
            ampl_y_mean = TGdatfile.get_moments("#04_ampl_factor_y", xcol=1, xs=3.0, xe=args.num_t_turb)['mean']
            ampl_z_mean = TGdatfile.get_moments("#05_ampl_factor_z", xcol=1, xs=3.0, xe=args.num_t_turb)['mean']
            ampl_factor = np.mean([ampl_x_mean, ampl_y_mean, ampl_z_mean])
            print("Adjusting time-averaged amplitudes ", [ampl_x_mean, ampl_y_mean, ampl_z_mean], color='cyan')
            print(" to mean amplitude ", ampl_factor, color='cyan')
            new_line = "ampl_factor        = "+str(ampl_factor)+"    ! Amplitude factor"
            cfp.replace_line_in_file(driv_gen_file, 'ampl_factor', new_line)
            new_line = "ampl_auto_adjust   = "+str(0)+"    ! Amplitude auto adjustment"
            cfp.replace_line_in_file(driv_gen_file, 'ampl_auto_adjust', new_line)

        if not args.adjust_amplitude:

            # edit flash.par
            flash_par_file = path+'/flash.par'
            new_line = "xmin                           = "+cfp.eform(-ic["L"]/2,6,print_leading_plus=True)
            cfp.replace_line_in_file(flash_par_file, 'xmin', new_line)
            new_line = "xmax                           = "+cfp.eform(+ic["L"]/2,6,print_leading_plus=True)
            cfp.replace_line_in_file(flash_par_file, 'xmax', new_line)
            new_line = "ymin                           = "+cfp.eform(-ic["L"]/2,6,print_leading_plus=True)
            cfp.replace_line_in_file(flash_par_file, 'ymin', new_line)
            new_line = "ymax                           = "+cfp.eform(+ic["L"]/2,6,print_leading_plus=True)
            cfp.replace_line_in_file(flash_par_file, 'ymax', new_line)
            new_line = "zmin                           = "+cfp.eform(-ic["L"]/2,6,print_leading_plus=True)
            cfp.replace_line_in_file(flash_par_file, 'zmin', new_line)
            new_line = "zmax                           = "+cfp.eform(+ic["L"]/2,6,print_leading_plus=True)
            cfp.replace_line_in_file(flash_par_file, 'zmax', new_line)
            new_line = "rho_ambient          =  "+cfp.eform(ic["rho"],6)
            cfp.replace_line_in_file(flash_par_file, 'rho_ambient', new_line)
            if star_formation_model:
                SFE_max = 0.3
                new_line = "sink_stop_on_total_mass   = "+cfp.eform(SFE_max*ic["mass"],6)
                cfp.replace_line_in_file(flash_par_file, 'sink_stop_on_total_mass', new_line)
                new_line = "sink_density_thresh       = "+cfp.eform(ic["sink_density"],6)
                cfp.replace_line_in_file(flash_par_file, 'sink_density_thresh', new_line)
                new_line = "sink_accretion_radius     = "+cfp.eform(ic["sink_radius"],6)
                cfp.replace_line_in_file(flash_par_file, 'sink_accretion_radius', new_line)
                new_line = "sink_softening_radius     = "+cfp.eform(ic["sink_radius"],6)
                cfp.replace_line_in_file(flash_par_file, 'sink_softening_radius', new_line)
            if star_formation_model:
                dumps_per_t_turb = 50
            else:
                dumps_per_t_turb = 10
            new_line = "tmax   = "+cfp.eform(ic["t_turb"]*args.num_t_turb,6)
            cfp.replace_line_in_file(flash_par_file, 'tmax', new_line)
            new_line = "dtmax  = "+cfp.eform(ic["t_turb"]/(2*dumps_per_t_turb),6)
            cfp.replace_line_in_file(flash_par_file, 'dtmax', new_line)
            new_line = "checkpointFileIntervalTime     = "+cfp.eform(ic["t_turb"],6)
            cfp.replace_line_in_file(flash_par_file, 'checkpointFileIntervalTime', new_line)
            new_line = "plotFileIntervalTime           = "+cfp.eform(ic["t_turb"]/dumps_per_t_turb,6)
            cfp.replace_line_in_file(flash_par_file, 'plotFileIntervalTime', new_line)
            new_line = "particleFileIntervalTime       = "+cfp.eform(ic["t_turb"]/dumps_per_t_turb,6)
            cfp.replace_line_in_file(flash_par_file, 'particleFileIntervalTime', new_line)
            new_line = "movie_dt_dump    = "+cfp.eform(ic["t_turb"]/dumps_per_t_turb,6)
            cfp.replace_line_in_file(flash_par_file, 'movie_dt_dump', new_line)
            if star_formation_model: # star formation (resim) model
                if m["resolution"] == 504:
                    iprocs = 6; jprocs = 8; kprocs = 8
                if m["resolution"] == 1024:
                    iprocs = 0; jprocs = 0; kprocs = 0
                if m["resolution"] == 4096:
                    iprocs = 0; jprocs = 0; kprocs = 0
            else: # driving model
                if m["resolution"] == 128:
                    iprocs = 2; jprocs = 2; kprocs = 2
                if m["resolution"] == 256:
                    iprocs = 4; jprocs = 4; kprocs = 4
                if m["resolution"] == 512:
                    iprocs = 8; jprocs = 8; kprocs = 8
                if m["resolution"] == 504:
                    iprocs = 3; jprocs = 4; kprocs = 4
                if m["resolution"] == 1008:
                    iprocs = 6; jprocs = 8; kprocs = 8
            new_line = "iProcs = "+str(iprocs)
            cfp.replace_line_in_file(flash_par_file, 'iProcs', new_line)
            new_line = "jProcs = "+str(jprocs)
            cfp.replace_line_in_file(flash_par_file, 'jProcs', new_line)
            new_line = "kProcs = "+str(kprocs)
            cfp.replace_line_in_file(flash_par_file, 'kProcs', new_line)
            # resim params
            if star_formation_model:
                new_line = "tinitial = "+cfp.eform(ic["t_turb"]*2,6)
                cfp.replace_line_in_file(flash_par_file, 'tinitial', new_line)
                new_line = "sim_length_scale = "+cfp.eform(ic["L"],6)
                cfp.replace_line_in_file(flash_par_file, 'sim_length_scale', new_line)
                new_line = "sim_dens_scale = "+cfp.eform(ic["rho"],6)
                cfp.replace_line_in_file(flash_par_file, 'sim_dens_scale', new_line)

            # link resim ics file
            if star_formation_model:
                path_to_driving_model = "../1008"+get_driving_str(m)+get_mach_str(m)+get_seed_str(m)
                cmd = "cd "+path+"; ln -s "+path_to_driving_model+"/OUT_hdf5_plt_cnt_0020 ./ics"
                print("now running '"+cmd+"'")
                sp_result = subprocess.run(cmd, shell=True)

        # edit job.sh
        job_file = path+'/job.sh'
        if not args.adjust_amplitude:
            new_line = "#PBS -N "+model_name
            cfp.replace_line_in_file(job_file, '#PBS -N', new_line)
            nprocs = iprocs*jprocs*kprocs
            if nprocs > 48: nprocs = int(np.ceil(nprocs/48))*48
            new_line = "#PBS -l ncpus="+str(nprocs)
            cfp.replace_line_in_file(job_file, '#PBS -l ncpus', new_line)
            new_line = "#PBS -l mem="+str(nprocs*4)+"GB"
            cfp.replace_line_in_file(job_file, '#PBS -l mem', new_line)
            line = cfp.find_line_in_file(job_file, '$PBS_NCPUS')
            new_line = line.replace('$PBS_NCPUS', str(iprocs*jprocs*kprocs))
            cfp.replace_line_in_file(job_file, '$PBS_NCPUS', new_line, search_str_position=None)
        else:
            line = cfp.find_line_in_file(job_file, 'shell.out')
            new_line = line.replace('shell.out_ampl_auto_adjust', 'shell.out01')
            cfp.replace_line_in_file(job_file, 'shell.out', new_line, search_str_position=None)

        if args.submit_job:
            cmd = "cd "+path+"; qsub job.sh"
            print("now running '"+cmd+"'", color='cyan')
            sp_result = subprocess.run(cmd, shell=True)
        """