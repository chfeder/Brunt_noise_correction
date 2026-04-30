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

# function to get simulation and dill filename(s)
def get_filename(res, void_size, dill=False):
    simpath = '../turb_foam_gamma_5_3_3d_'+str(res)+'_'+str(void_size)+'_bubbles/'
    filename = simpath+'NIF_hdf5_plt_cnt_0375'
    if dill: filename += "_turb_data.dill"
    return filename

# figure output path
figpath = '../figs/'
