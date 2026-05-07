import cfpack as cfp
from cfpack import stop, print
cfp.load_plot_style()
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import h5py
from collections import defaultdict
from types import SimpleNamespace
import dill

# === global variables and functions ===

# centre, left, right, back, front box
BOX_list = ['c_box', 'l_box', 'r_box', 'b_box', 'f_box']

# 3 different projections (LOS)
LOS_list = [0, 1, 2] # x, y, z

# SNR list (8 values, with the 0th index holding near infinite SNR)
SNR_list = [1e6, 0.01, 0.1, 1, 3, 5, 10, 100]

# global void size strings and formats
void_sizes = ['12p5mu', '50mu', '100mu']
void_colors = {'12p5mu':'magenta', '50mu':'blue', '100mu':'green'}
void_markers = {'12p5mu':'s', '50mu':'o', '100mu':'D'}
void_size_labels = {'12p5mu':r'$12.5\,\mu\mathrm{m}$',
                    '50mu':r'$50\,\mu\mathrm{m}$',
                    '100mu':r'$100\,\mu\mathrm{m}$'}

# function to get simulation and dill filename(s)
def get_filename(res, void_size, noise_exp=1.0, res_resample=1.0, dill=False):
    simpath = '../turb_foam_gamma_5_3_3d_'+str(res)+'_'+str(void_size)+'_bubbles/'
    filename = simpath+'NIF_hdf5_plt_cnt_0375'
    if dill:
        noise_str = ""
        if noise_exp != 1.0:
            noise_str = "_noise_exp_"+str(noise_exp)
        res_resample_str = ""
        if res_resample != 1.0:
            res_resample_str = "_res_resample_"+str(res_resample)
        filename += "_turb_data"+noise_str+res_resample_str+".dill"
    return filename

# figure output path
figpath = '../figs/'
