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

# 3 different projections (LOS)
LOS_list = [0, 1, 2] # x, y, z

# figure output path
figpath = '../figs/'

# ---------------- 5 BOX SIMULATION ONLY ----------------

# centre, left, right, back, front box
BOX_list = ['c_box', 'l_box', 'r_box', 'b_box', 'f_box']

# SNR list (8 values, with the 0th index holding near infinite SNR)
SNR_list = [1e6, 0.01, 0.1, 1, 3, 5, 10, 100]

# global void size strings and formats
void_sizes = ['12p5mu', '50mu', '100mu']
void_colors = {'12p5mu':'magenta', '50mu':'blue', '100mu':'green'}
void_linestyles = {'12p5mu':'dotted', '50mu':'solid', '100mu':'dashed'}
void_markers = {'12p5mu':'s', '50mu':'o', '100mu':'D'}
void_size_labels = {'12p5mu':r'$12.5\,\mu\mathrm{m}$',
                    '50mu':r'$50\,\mu\mathrm{m}$',
                    '100mu':r'$100\,\mu\mathrm{m}$'}

# global resolution resampling (k_max) variables
res_resample_list = [0.25, 0.5, 1.0, 2.0, 4.0]
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list('magenta_scale', ['#f7d6f7', '#ff00ff', '#4b004b'])
res_colors_magenta = cmap(np.linspace(0, 1, 5))
res_colors_blue = plt.cm.Blues(np.linspace(0.3, 1.0, 5))
res_colors_green = plt.cm.Greens(np.linspace(0.3, 1.0, 5))

# ------------------ UNIFORM SIMULATION ONLY ----------------

# Global settings
resolutions = [256]                                         # Simulation resolution
flow_types = ['c', 0.1, 0.2, 0.3, 0.4, 'm', 0.75, 's']      # s=Soluniodal (zeta=1), m=Mixed (zeta=0.5), c=Compressive (zeta=0)
flow_types_nums = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1])
mach_nums = [0.1, 0.2, 0.5, 1, 2, 5, 10]                    # Mach numbers
project_path = "/scratch/ek9/ag9727/asc/"                   # Path to folder
min_time = 30                                               # Files correponding to times before this will be discarded
max_time = 100                                              # Likewise, after this will be discarded
SNR_list_uniform = [1e6, 0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 100000] # Signal to Noise Ratios

# Test case:
"""
resolutions = [256]                         # simulation resolution (default first)
flow_types = ['s']                          # Soluniodal, Mixed, Compressive
mach_nums = [0.1]                           # Mach numbers
project_path = "/scratch/ek9/ag9727/asc/"   # Path to folder
min_time = 30                               # Files correponding to times before this will be discarded
max_time = 31                               # Likewise, after this will be discarded
"""


flow_colors = {
    's':'magenta', 'c':'blue', 'm':'green',
    0.1:'magenta', 1:'blue', 10:'green',
}
flow_linestyles = {
    's':'dotted', 'c':'dotted', 'm':'dotted',
    0.1:'dotted', 1:'dotted', 10:'dotted',
}
flow_markers = {
    's':'o', 'c':'o', 'm':'o',
    0.1:'o', 1:'o', 10:'o',
}
flow_labels = {
    's':'Solenoidal', 'c':'Compressive', 'm':'Mixed',
    0.1: 0.1, 1:1, 10:10,
}

flow_to_num = dict()
for x in flow_types_nums:
    flow_to_num[x] = x
flow_to_num['s'] = 1
flow_to_num['m'] = 0.5
flow_to_num['c'] = 0

zeta_colours = plt.cm.magma(np.linspace(0, 0.9, len(flow_types)))

mach_colors = plt.cm.viridis(np.linspace(0, 0.9, len(mach_nums)))
mach_linestyles = {0.1:'dotted', 0.2:'dotted', 0.5:'dotted', 1:'dotted', 10:'dotted'}
mach_markers = {0.1:'o', 0.2:'o', 0.5:'o', 1:'o', 10:'o'}
mach_labels = {
    0.1: r'$\mathcal{M}=0.1$', 1:r'$\mathcal{M}=1$', 10:r'$\mathcal{M}=10$',
}


# ========== global functions ============

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

# ===================================================================
# Given a side length N and noise_exp, return the noise values over an N x N grid
# The noise values are normalised to have a mean of 0 and std of 1
def get_noise(N, noise_exp=1.0, kmin=None, kmax=None, std=1.0, random_phases=True, random_amplitudes=True, random_seed=None):
    # noise_exp = 1 creates a white-noise spectrum. The noise_exp is called 'alpha' in the paper (cf. Sec. 4.1).
    # Note that alpha = beta + 1, and hence, beta = alpha - 1. With alpha = 1, we have beta = 0, which is white noise.
    # The 1D azimuthally-integrated spectrum scales as k^noise_exp for 2D fields, and since we plot it as such,
    # we use alpha=noise_exp to denote the noise exponent in the paper.
    isovec = [1,1,1]
    kI = np.fft.fftfreq(N)
    kx, ky = np.meshgrid(kI, kI, indexing='ij') # 2D
    rrr = np.sqrt(isovec[0]*kx**2+isovec[1]*ky**2) # 2D
    Ahat = np.zeros_like(rrr*1j)
    k = np.unique(np.abs(kx))
    if kmin is None:
        kmin = 0
    if kmax is None:
        kmax = N//2
    rmin = k[kmin]
    rmax = k[kmax]
    ok = (rrr>rmin) & (rrr<=rmax)

    beta = noise_exp - 1.0 # beta in Eq. (10) in Sec. 4.1
    betaprime = beta / 2
    Ahat[ok] = rrr[ok]**betaprime
    if random_seed:
        print('Setting random seed ', random_seed)
        np.random.seed(random_seed)
    if random_amplitudes:
        rn = np.random.random(Ahat.size)
        rn.shape = Ahat.shape
        Ahat *= rn
    Ahat[0,0] = 0 # mean amplitude, 2D
    if random_phases:
        phi = (np.random.random(Ahat.size)-0.5)*np.pi*2
        phi.shape = Ahat.shape
        Ahatmag = np.abs(Ahat)
        Ahat = Ahatmag*np.cos(phi)+Ahatmag*np.sin(phi)*1j
        Ahat[0,0] = np.abs(Ahat[0,0]) # 2D
    noise = np.fft.irfftn(Ahat, s=Ahat.shape)
    noise -= noise.mean() # correct tiny deviations from 0 mean
    noise /= noise.std() # normalise to target STD
    noise *= std # normalise to target STD
    return noise

# ===================================================================
# Given a signal (N x N array), add noise to it according to the desired SNR
def add_noise(signal, SNR, noise_exp=1.0):
    N = signal.shape[0]
    noise = get_noise(N, noise_exp) # return noise with 0 mean and std of 1
    scale = np.std(signal) / SNR
    noise *= scale
    noisy_data = signal + noise
    return noisy_data, noise

# ===================================================================
def get_Brunt_data(k, P, k_noise=np.inf): # k_noise achieves denoising
    ind = (k > 0) & (k <= k_noise)
    P2D = P[ind].sum()
    P3D = (2*k[ind]*P[ind]).sum()
    Brunt_R_to_half = np.sqrt(P2D/P3D)
    Brunt_sigma_2D = np.sqrt(P2D/P[0]) # normalised by the mean (np.sqrt(P[0]))
    Brunt_sigma_3D = Brunt_sigma_2D / Brunt_R_to_half
    return Brunt_sigma_3D, Brunt_R_to_half

# ===================================================================
# Given a signal (N x N array), and an SNR, add noise to it for a specified number of trials
# Return the resulting data in the form of a ret class, defined below
def noise_iterations(iterations, signal, SNR, noise_exp=1.0):
    # prepare output lists
    sigma_2D_list = []
    PS_SpN_list = [] # PS = Power Spectrum, SpN = Signal plus Noise
    PS_N_list = []
    k_noise_list = []
    PS_S_at_k_noise_list = []
    Brunt_sigma_3D_list = []
    Brunt_R_to_half_list = []
    Brunt_sigma_3D_denoised_list = []
    Brunt_R_to_half_denoised_list = []
    # signal spectrum
    spect_S = cfp.get_spectrum(signal, verbose=0)
    k = spect_S['k']
    PS_S = spect_S['P_tot']
    # loop over noise iterations
    for iteration in range(iterations):
        # add noise
        noisy_data, noise = add_noise(signal, SNR, noise_exp)
        # noisy sigma_2D
        sigma_2D = (noisy_data/noisy_data.mean()).std()
        sigma_2D_list.append(sigma_2D)
        # compute spectra
        PS_SpN = cfp.get_spectrum(noisy_data, verbose=0)['P_tot']
        PS_SpN_list.append(PS_SpN)
        PS_N = cfp.get_spectrum(noise, verbose=0)['P_tot']
        PS_N_list.append(PS_N)
        # k_noise and PS_S value at k_noise


        if not all(PS_S[1:] > PS_N[1:]):
            # Just this on its own (without the if) can cause k_noise to
            # be set to a local min should there be an uptick in the far right
            k_noise_ind = np.argmin(np.abs(np.log10(PS_S[1:]/PS_N[1:])))
        else:
            k_noise_ind = len(k[1:]) - 1

        k_noise = k[1:][k_noise_ind]
        k_noise_list.append(k_noise)
        PS_S_at_k_noise_list.append(PS_S[1:][k_noise_ind])
        # Brunt sigma_3D and R^(1/2) without noise correction
        Brunt_sigma_3D, Brunt_R_to_half = get_Brunt_data(k, PS_SpN, k_noise=np.inf)
        Brunt_sigma_3D_list.append(Brunt_sigma_3D)
        Brunt_R_to_half_list.append(Brunt_R_to_half)
        # Brunt sigma_3D and R^(1/2) denoised
        Brunt_sigma_3D_denoised, Brunt_R_to_half_denoised = get_Brunt_data(k, PS_SpN, k_noise)
        Brunt_sigma_3D_denoised_list.append(Brunt_sigma_3D_denoised)
        Brunt_R_to_half_denoised_list.append(Brunt_R_to_half_denoised)
    class ret:
        def __init__(self, SNR_, noise_exp_, signal_, noise_, noisy_data_, sigma_2D_,
                        k_, PS_S_, PS_N_, PS_SpN_, k_noise_, PS_S_at_k_noise_,
                        Brunt_sigma_3D_, Brunt_R_to_half_,
                        Brunt_sigma_3D_denoised_, Brunt_R_to_half_denoised_):
            self.SNR = np.array([SNR_]) # add extra ghost dim to allow simple stats function below
            self.noise_exp = np.array([noise_exp_])
            self.signal = np.array([signal_])
            self.noise = np.array([noise_])
            self.noisy_data = np.array([noisy_data_])
            self.sigma_2D = np.array(sigma_2D_)
            self.k = np.array([k_])
            self.PS_S = np.array([PS_S_]) # add extra ghost dim to allow simple stats function below
            self.PS_N = np.array(PS_N_)
            self.PS_SpN = np.array(PS_SpN_)
            self.k_noise = np.array(k_noise_)
            self.PS_S_at_k_noise = np.array(PS_S_at_k_noise_)
            self.Brunt_sigma_3D = np.array(Brunt_sigma_3D_)
            self.R_to_half = np.array(Brunt_R_to_half_)
            self.Brunt_sigma_3D_denoised = np.array(Brunt_sigma_3D_denoised_)
            self.R_to_half_denoised = np.array(Brunt_R_to_half_denoised_)
    return ret(SNR, noise_exp, signal, noise, noisy_data, sigma_2D_list,
                k, PS_S, PS_N_list, PS_SpN_list, k_noise_list, PS_S_at_k_noise_list,
                Brunt_sigma_3D_list, Brunt_R_to_half_list,
                Brunt_sigma_3D_denoised_list, Brunt_R_to_half_denoised_list)

# ===================================================================
# --- function that computes statistics over the noise samples ---
# Given a list of ret objects (where the list is all ret objects for a particular SNR,
# and each ret corresponds to a single projection) return the list [avg, std, p16, p84],
# where each entry is a ret object, containing all the data for that statiscal value.
def get_stats_object(objs):
    # -----------------------------------------------------------
    # Given an array containing all the data for a single attribute (across all ret objects),
    # return the array [avg, std, 16th, 84th] for it. Some may not be scalars.
    def get_stats(arr, log_stats=True):
        if log_stats:
            arr[arr==0] = 1e-99 # fix 0 values in log statistics
            data_arr = np.log10(arr) # log statistics
        else:
            data_arr = arr # linear statistics
        stats = [np.nan, np.nan, np.nan, np.nan] # avg, std, 16th, 84th
        # statistics over axes 0 and 1 = boxes and noise realisations
        stats[0] = np.mean      (data_arr,     axis=(0, 1)) # avg
        stats[1] = np.std       (data_arr,     axis=(0, 1)) # std
        stats[2] = np.percentile(data_arr, 16, axis=(0, 1)) # 16th
        stats[3] = np.percentile(data_arr, 84, axis=(0, 1)) # 84th
        if log_stats: # log statistics
            for i in range(len(stats)):
                stats[i] = 10**stats[i]
        return stats
    # -----------------------------------------------------------

    # Make a list of 4 new ret objects to return, the avg, std, p16, p84
    cls = type(objs[0])
    ret_stats_data = [cls.__new__(cls) for _ in range(4)]  # avg, std, p16, p84

    # For each (attribute name, attribute data) in the ret object...
    for k, v in vars(objs[0]).items():
        # print('Computing statistics on '+k)

        # Get every value for the single attribute k for every ret object in a single list
        vals = [getattr(o, k) for o in objs]
        arr = np.asarray(vals)

        # Get all the statistics (the list [avg, std, p16, p84]) for the attribute
        if k in ['noise_exp', 'noise', 'noisy_data', 'k']:
            log_stats = False
        else:
            log_stats = True
        stats = get_stats(arr, log_stats)

        # Store all the results into the ret_stats_data
        for i in range(len(stats)):
            setattr(ret_stats_data[i], k, stats[i])

    return ret_stats_data

# ===================================================================
def make_mach_str(mach: float) -> str:
    if mach == int(mach):
        mach_str = 'M' + str(int(mach))
    else:
        mach_str = 'M' + str(mach).replace('.', 'p')
    return mach_str

def make_flow_str(flow: str | float | int) -> str:
    if flow in ['s', 'm', 'c']:
        return flow
    else:
        # zeta (which is in [0,1]) refers to a blend of sol and comp
        zeta = float(flow)
        if zeta == int(zeta):
            return 'z' + str(int(zeta))
        else:
            return 'z' + str(zeta).replace('.', 'p')

# ---- generates the model as a string for a given simulation, ie, "N256cM0p2", given the below params
def make_model_name_params(res, flow_type: str | float | int, mach) -> str:
    res = str(res)

    flow_str = make_flow_str(flow_type)

    mach = float(mach)
    mach_str = make_mach_str(mach)

    return 'N' + res + flow_str + mach_str

