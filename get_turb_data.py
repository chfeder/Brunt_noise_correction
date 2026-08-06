#!/usr/bin/env python

from globals_turb_data import *
from FlashFile2d import FlashFile2d
import flashlib

# ===================================================================
def get_turb_box(filename, box, res_resample=1.0):
    return flashlib.FlashGG(filename).GetUniformGrid()

    """
    import yt
    yt.funcs.mylog.setLevel("ERROR")
    # box dimensions
    delta_a = 0.02
    edge_mtop = 0.35+delta_a
    edge_mbot = 0.35-delta_a
    # read FLASH file and extract box region via yt
    ds = yt.load(filename)
    if box == 'c_box':
       ledge = ds.arr([-delta_a, edge_mbot, -delta_a],'code_length')
       redge = ds.arr([+delta_a, edge_mtop, +delta_a],'code_length')
    if box == 'r_box':
        ledge = ds.arr([delta_a, edge_mbot, -delta_a],'code_length')
        redge = ds.arr([delta_a+2*delta_a, edge_mtop, delta_a],'code_length')
    if box == 'f_box':
        ledge = ds.arr([-delta_a, edge_mbot, delta_a],'code_length')
        redge = ds.arr([+delta_a, edge_mtop, delta_a+2*delta_a],'code_length')
    if box == 'l_box':
        ledge = ds.arr([-delta_a-2*delta_a, edge_mbot, -delta_a],'code_length')
        redge = ds.arr([-delta_a, edge_mtop, delta_a],'code_length')
    if box == 'b_box':
        ledge = ds.arr([-delta_a, edge_mbot, -delta_a-2*delta_a],'code_length')
        redge = ds.arr([+delta_a, edge_mtop, -delta_a],'code_length')
    dx = ds.index.get_smallest_dx()
    size = ((redge-ledge)/dx).astype('int')
    level = 0
    section = ds.covering_grid(level, ledge, size)
    dens = section['density'].value
    # requested resolution resampling by res_resample factor
    if res_resample != 1.0:
        new_shape = tuple(np.round(np.array(dens.shape)*res_resample).astype(int))
        print('Requested -- reshaping ', dens.shape, ' -> ', new_shape, warn=True)
        mean, std = dens.mean(), dens.std()
        dens = cfp.congrid(dens, new_shape)
        dens = (dens - dens.mean())/dens.std()*std + mean # correct to get exact previous mean and std
    # treat odd shape by interpolation
    if (np.array(dens.shape)%2)[0] != 0:
        new_shape = tuple(np.array(dens.shape)-1)
        print('Odd shape -- reshaping ', dens.shape, ' -> ', new_shape, warn=True)
        mean, std = dens.mean(), dens.std()
        dens = cfp.congrid(dens, new_shape)
        dens = (dens - dens.mean())/dens.std()*std + mean # correct to get exact previous mean and std
    return dens # density on uniform grid
    """

# ===================================================================
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
    if kmin is None: kmin = 0
    if kmax is None: kmax = N//2
    rmin = k[kmin]
    rmax = k[kmax]
    ok = (rrr>rmin)*(rrr<=rmax)
    ok_r = ok
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
def noise_iterations(iterations, signal, SNR, noise_exp=1.0):
    # prepare output lists
    sigma_2D_list = []
    PS_SpN_list = []
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
        k_noise_ind = np.argmin(np.abs(np.log10(PS_S[1:]/PS_N[1:])))
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
def get_stats_object(objs):
    # -----------------------------------------------------------
    def get_stats(arr, log_stats=True):
        if log_stats:
            arr[arr==0] = 1e-99 # fix 0 values in log statistics
            data_arr = np.log10(arr) # log statistics
        else: data_arr = arr # linear statistics
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
    cls = type(objs[0])
    ret_stats_data = [cls.__new__(cls) for _ in range(4)]  # avg, std, p16, p84
    for k, v in vars(objs[0]).items(): # loop over all keys
        # print('Computing statistics on '+k)
        vals = [getattr(o, k) for o in objs]
        arr = np.asarray(vals)
        if k in ['noise_exp', 'noise', 'noisy_data', 'k']: log_stats = False
        else: log_stats = True
        stats = get_stats(arr, log_stats)
        for i in range(len(stats)):
            setattr(ret_stats_data[i], k, stats[i])
    return ret_stats_data
# -----------------------------------------------------------


# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    # user changes here
    resolutions = [1024, 512, 2048] # simulation resolution (default first)
    void_sizes = ['50mu', '12p5mu', '100mu'] # foam void sizes (default first)
    noise_samples = 100 # convergence was checked also with 1000 noise samples, which yields practically the same
    print('Number of noise samples: '+str(noise_samples), color='lightred')

    timer = cfp.timer("MAIN")

    # loop over simulation resolutions
    for res in resolutions:

        # loop over foam void sizes
        for void_size in void_sizes:
            if res != 1024 and void_size != '50mu': continue # skip 512 and 2048 resolutions for 12.5 and 100 mu
            print("===== SIMULATION RESOLUTION AND VOID SIZE = "+str(res)+", "+void_size+" =====", color='lightred')

            # loop over noise exponents
            for noise_exp in [1.0, 0.0, 2.0, 0.5, 1.5]:

                # skip noise exponents other than 1 (white noise) for foam sizes 12.5mu and 100mu
                if noise_exp != 1.0 and (void_size != '50mu' or res != 1024): continue

                noise_exp_info = str(noise_exp)
                if noise_exp == 1.0: noise_exp_info += ' (white noise)'
                print("===== NOISE EXPONENT = "+noise_exp_info+" =====", color='cyan')

                # loop over resolution resampling cases
                for res_resample in res_resample_list:
                    if res_resample != 1.0 and res != 1024: continue

                    if res_resample != 1.0:
                        print("===== RESOLUTION RESAMPLING FACTOR = "+str(res_resample)+" =====", color='yellow')

                    # data collection per box and/or LOS
                    sigma_3D = []
                    data = [[] for _ in range(len(SNR_list))]

                    filename = get_filename(res, void_size)

                    # loop over the 5 boxes
                    for box in BOX_list:

                        dens = get_turb_box(filename, box, res_resample=res_resample)
                        print("=== Read "+box+" with dims "+str(dens.shape)+" from file "+filename, color='lightblue')

                        sigma_3D = (dens/dens.mean()).std()

                        # loop over the 3 projection axes (x,y,z)
                        for LOS in LOS_list:
                            print("computing projection along axis "+str(LOS), color='lightblue')
                            # five box locations
                            delta_a = 0.02
                            proj = dens.sum(axis=LOS) / dens.shape[LOS] * 2*delta_a

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
                    dillfile = get_filename(res, void_size, noise_exp=noise_exp,
                                            res_resample=res_resample, dill=True)
                    with open(dillfile, "wb") as f:
                        dill.dump(turb_data, f)
                    print("Wrote '"+dillfile+"'.", color='magenta')

    timer.report()