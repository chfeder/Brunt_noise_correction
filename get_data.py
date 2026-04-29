#!/usr/bin/env python

from imports import *

# === user changes here ===
res = '1024' # resolution
void_sizes = ['12p5mu', '50mu', '100mu']
void_sizes = ['50mu']
noise_samples = 2
print('Number of noise samples: '+str(noise_samples), color='lightred')

# box dimensions
delta_a = 0.02
edge_mtop = 0.35+delta_a # 0.3675
edge_mbot = 0.35-delta_a # 0.3275


# ===================================================================
def get_box(filename, box):
    import yt
    yt.funcs.mylog.setLevel("ERROR")
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
    # treat odd shape by interpolation
    if (np.array(dens.shape)%2)[0] != 0:
        new_shape = tuple(np.array(dens.shape)-1)
        print('odd shape -- reshaping from -> to: ', dens.shape, ' -> ', new_shape, warn=True)
        total = dens.sum()
        dens = cfp.congrid(dens, new_shape)
        dens *= total/dens.sum()
    return dens # density on uniform grid

# ===================================================================
def get_noise(N, noise_exp=0.0, kmin=None, kmax=None, std=1.0, random_phases=True, random_amplitudes=True, random_seed=None):
    # noise_exp = 0.0 creates a white-noise spectrum.
    # The 1D azimuthally-integrated spectrum scales as k^(noise_exp+1) for 2D fields.
    isovec=[1,1,1]
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
    betaprime = noise_exp / 2
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
    noise /= noise.std() # normalise to target STD
    noise *= std # normalise to target STD
    return noise

# ===================================================================
def add_noise(signal, SNR, noise_exp=0.0):
    N = signal.shape[0]
    noise = get_noise(N, noise_exp)
    scale = np.std(signal) / SNR
    scale /= noise.std()
    noise *= scale
    noisy_data = signal + noise
    return noisy_data, noise

# ===================================================================
def get_Brunt_data(k, P, k_noise):
    ind = (k > 0) & (k <= k_noise)
    P2D = P[ind].sum()
    P3D = (2*k[ind]*P[ind]).sum()
    Brunt_R_to_half = np.sqrt(P2D/P3D)
    Brunt_sigma_2D = np.sqrt(P2D/P[0]) # normalised by the mean (np.sqrt(P[0]))
    Brunt_sigma_3D = Brunt_sigma_2D / Brunt_R_to_half
    return Brunt_sigma_3D, Brunt_R_to_half

# ===================================================================
def noise_iterations(iterations, signal, SNR, noise_exp=0.0):
    # prepare output lists
    sigma_2D_list = []
    PS_SpN_list = []
    PS_N_list = []
    k_noise_list = []
    Brunt_sigma_3D_list = []
    Brunt_R_to_half_list = []
    # signal spectrum
    spect_S = cfp.get_spectrum(signal, verbose=0)
    k = spect_S['k']
    PS_S = spect_S['P_tot']
    # loop over noise iterations
    for iteration in range(iterations):
        # add noise
        noisy_data, noise = add_noise(signal, SNR)
        # noisy sigma_2D
        sigma_2D = (noisy_data/noisy_data.mean()).std()
        sigma_2D_list.append(sigma_2D)
        # compute spectra
        PS_SpN = cfp.get_spectrum(noisy_data, verbose=0)['P_tot']
        PS_SpN_list.append(PS_SpN)
        PS_N = cfp.get_spectrum(noise, verbose=0)['P_tot']
        PS_N_list.append(PS_N)
        # k_noise
        k_noise_ind = np.argmin(np.abs(np.log10(PS_S[1:]/PS_N[1:])))
        k_noise = k[1:][k_noise_ind]
        k_noise_list.append(k_noise)
        # Brunt sigma_3D and R^(1/2)
        Brunt_sigma_3D, Brunt_R_to_half = get_Brunt_data(k, PS_SpN, k_noise)
        Brunt_sigma_3D_list.append(Brunt_sigma_3D)
        Brunt_R_to_half_list.append(Brunt_R_to_half)
    class ret:
        def __init__(self, SNR_, noise_exp_, signal_, noise_, noisy_data_,
                        k_, PS_S_, PS_N_, PS_SpN_, k_noise_,
                        Brunt_sigma_3D_, Brunt_R_to_half_, sigma_2D_):
            self.SNR = np.array([SNR_]) # add extra ghost dim to allow simple stats function below
            self.noise_exp = np.array([noise_exp_])
            self.signal = np.array([signal_])
            self.noise = np.array([noise_])
            self.noisy_data = np.array([noisy_data_])
            self.k = np.array([k_]) 
            self.PS_S = np.array([PS_S_]) # add extra ghost dim to allow simple stats function below
            self.PS_N = np.array(PS_N_)
            self.PS_SpN = np.array(PS_SpN_)
            self.k_noise = np.array(k_noise_)
            self.Brunt_sigma_3D = np.array(Brunt_sigma_3D_)
            self.R_to_half = np.array(Brunt_R_to_half_)
            self.sigma_2D = np.array(sigma_2D_)
    return ret(SNR, noise_exp, signal, noise, noisy_data,
                k, PS_S, PS_N_list, PS_SpN_list, k_noise_list,
                Brunt_sigma_3D_list, Brunt_R_to_half_list, sigma_2D_list)


# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    timer = cfp.timer("MAIN")

    for void_size in void_sizes:

        print("===== NOW WORKING ON VOID SIZE "+void_size+" =====", color='lightred')

        # BOX OPTIONS: 15 boxes per void size.
        # 3 different box locations for the x LOS: left, center, right. and 2 more for the z LOS: front & back
        boxes = ['c_box', 'l_box', 'r_box', 'b_box', 'f_box']

        # 3 different projections (LOS)
        LOS_list = [0, 1, 2] # x, y, z

        SNR_list = [1e6, 0.01, 0.1, 1, 3, 5, 10, 100]

        # data collection per box and/or LOS
        sigma_3D = []
        data = [[] for _ in range(len(SNR_list))]

        # 5 boxes
        for box in boxes:

            filename = '../turb_foam_gamma_5_3_3d_'+str(res)+'_'+str(void_size)+'_bubbles/NIF_hdf5_plt_cnt_0375'
            print("=== Reading "+box+" from file "+filename, color='lightblue')
            dens = get_box(filename, box)
            sigma_3D = (dens/dens.mean()).std()

            # 3 projection axes
            for LOS in LOS_list:
                print("computing projection along axis "+str(LOS), color='lightblue')
                # five box locations
                proj = dens.sum(axis=LOS) / dens.shape[LOS] * 2*delta_a

                # 7 SNR
                for iSNR, SNR in enumerate(SNR_list):

                    ret = noise_iterations(noise_samples, proj, SNR)
                    ret.sigma_3D = np.array([sigma_3D]) # attach sigma_3D
                    data[iSNR].append(ret)

                    if 0:

                        if target != 1000000:
                            rsqrt_k = []
                            medspec2D = collectspec
                            medspec2D_ave = collectspec_ave
                            medspec3D = 2 * kvals2D * medspec2D
                            medspec3D_ave = 2 * kvals2D * medspec2D_ave
                            # interpolation to get accurate intersection k values
                            signalonly = collect_spectras[1000000][-1]
                            meannoisespec = collectnoisespec
                            #kcross, P_at_kcross = get_k_cross(kvals, meannoisespec, signalonly)
                            # get intersection
                            # print('Target SNR = ', target, color='yellow')
                            #kcross = np.argmin(np.abs(signalonly_new - meannoisespec_new))
                            #collect_p2dsInt[target].append(signalonly_new[kcross])
                            #kcross *= kvals2D[-1] / 500 # scale back to original k grid
                            # print ('S-N intersection k value = ', kcross, color='yellow')
                            #if kcross > 1: stop()

                            #if 1:
                            #    _, sigmaNdenoised, rsqrt_true, sigmadenoised = parsevals.parsevals_th(1, medspec2D, medspec3D, \
                            #        kvals2D, sigmaN_sq=sigmaN_med, the_k=kcross, p2d_ave=medspec2D_ave, p3d_ave=medspec3D_ave)
                            #rsqrt_ratio = rsqrt_true / rsqrt_med
                            kcross = 66
                            rsqrt_true = 66
                            rsqrt_ratio = 66
                            sigmaNdenoised = 66
                            sigmadenoised = 66

                            collect_sigmasNdenoised[target].append(sigmaNdenoised)
                            collect_sigmasdenoised[target].append(sigmadenoised)
                            collect_kvals[target].append(kcross)
                            collect_rsqrttrue[target].append(rsqrt_true)
                            collect_rsqrtsratio[target].append(rsqrt_ratio)
                            collect_sigmas[target].append(sigma_med)
                            collect_rsqrt[target].append(rsqrt_med)

        # -----------------------------------------------------------
        def stats_object(objs):
            # -----------------------------------------------------------
            def get_stats(arr, log_stats=True):
                if log_stats: data_arr = np.log10(arr) # log statistics
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
            ret_stats_data = [cls.__new__(cls)]*4 # avg, std, 16th, 84th
            for k, v in vars(objs[0]).items(): # loop over all keys
                print('Computing statistics on '+k)
                vals = [getattr(o, k) for o in objs]
                arr = np.asarray(vals)
                if k in ['noise_exp', 'noise', 'k']: log_stats = False
                else: log_stats = True
                stats = get_stats(arr, log_stats)
                for i in range(len(stats)):
                    setattr(ret_stats_data[i], k, stats[i])
            return ret_stats_data
        # -----------------------------------------------------------

        avg, std, p16, p84 = stats_object(data[0])

    # plot example maps
    for isnr in [2,3,4,6]:
        # check SNR match
        LOS = 1 # y axis
        dat = data[isnr][LOS]
        if dat.SNR[0] != SNR_list[isnr]: stop()
        cfp.plot_map(dat.noisy_data[0], axes_format=["", ""], vmin=5e-3, vmax=1e-2, show=True)
        stop()


        # convert the default dict to a vector array
        if not no_noise:
            h5file = '../h5files/15cubes15spectra_i1k_1p0_'+void_size+'_'+str(res)+'.h5'
            Fptr = h5py.File(h5file, 'w')
            Fptr['sigmas'] = [*collect_sigmas.values()]
            Fptr['rsqrtfull'] = [*collect_rsqrt.values()]
            Fptr['kvals'] = [*collect_kvals.values()]
            Fptr['rsqrttrue'] = [*collect_rsqrttrue.values()]
            Fptr['ratios'] = [*collect_rsqrtsratio.values()]
            Fptr['kvals2D'] = kvals2D

            Fptr['powintersects'] = [*collect_p2dsInt.values()]
            Fptr['columnddenoise'] = [*collect_sigmasNdenoised.values()]
            Fptr['sigmasdenoised'] = [*collect_sigmasdenoised.values()]

            Fptr['thepowspec'] = [*collect_spectras.values()]
            Fptr['thepowspec_ave'] = [*collect_spectras_ave.values()]
            Fptr['thenoisepowspec'] = [*collect_noisespectras.values()]
            Fptr.close()
            print(h5file+' written.', color='magenta')

        if 0:
            Fptr['sigma_3Ds'] = sigmas3D

        if no_noise:
            h5file = '../h5files/15nonoise_%s.h5'%void_size
            Fptr = h5py.File(h5file, 'w')
            Fptr['sigma_brunts'] = sigmaBrunt
            Fptr['sigma_2Ds'] = sigma2D
            Fptr['sigma_rhalfs'] = rhalfs
            Fptr.close()
            print(h5file+' written.', color='magenta')





        # IMAGES

        if 0:

            # LOAD raw DATA
            ds = yt.load('../turb_foam_gamma_5_3_3d_'+str(res)+'_'+str(void_size)+'_bubbles/NIF_hdf5_plt_cnt_0375')

            proj_2 = yt.SlicePlot(ds, 2, ('gas', 'density'))
            proj_2.annotate_line((-3*delta_a, edge_mtop, 0), (3*delta_a, edge_mtop, 0), coord_system='data')
            proj_2.annotate_line((-3*delta_a, edge_mbot, 0), (3*delta_a, edge_mbot, 0), coord_system='data')
            proj_2.annotate_line((-3*delta_a, edge_mtop, 0), (-3*delta_a, edge_mbot, 0), coord_system='data')
            proj_2.annotate_line((3*delta_a, edge_mtop, 0), (3*delta_a, edge_mbot, 0), coord_system='data')
            proj_2.annotate_line((-delta_a, edge_mtop, 0), (-delta_a, edge_mbot, 0), coord_system='data')
            proj_2.annotate_line((delta_a, edge_mtop, 0), (delta_a, edge_mbot, 0), coord_system='data')
            proj_2.set_cmap(('gas','density'), cmap='magma')
            colorbar_2 = proj_2.plots['density'].cb
            vmin = 1e-3
            vmax = 1e1
            colorbar_2.vmin = vmin
            colorbar_2.vmax = vmax
            proj_2.set_axes_unit('mm')
            proj_2.set_font_size(20)
            proj_2.set_ylabel(r'$y \,\mathrm{(mm)}$')
            proj_2.set_xlabel('')
            proj_2.plots['gas', 'density'].axes.get_xaxis().set_tick_params(labelbottom=False)
            proj_2.set_zlim('density', vmin, vmax)
            proj_2.hide_colorbar()
            proj_2.save('z_slice.pdf')

            dx = 3.9062486030E-04
            projection_length = 2*dx  # try with 1dx
            ledge = ds.arr([-0.15, edge_mbot, -0.15],'code_length')
            ledge_slice = ds.arr([-0.15, edge_mcenter-projection_length/2, -0.15],'code_length')
            redge = ds.arr([0.15, edge_mtop, 0.15],'code_length')
            redge_slice = ds.arr([0.15, edge_mcenter+projection_length/2, 0.15],'code_length')
            center =0.5*ledge + 0.5*redge
            the_region_slice = ds.region(center, ledge_slice, redge_slice)  
            the_region = ds.region(center, ledge, redge)
            proj_1 = yt.ProjectionPlot(ds, 1, ('gas','density'), data_source=the_region_slice) # slices don't accept CGS
            proj_1.annotate_line((-delta_a, 0, -3*delta_a), (-delta_a, 0, 3*delta_a), coord_system='data')
            proj_1.annotate_line((delta_a, 0, -3*delta_a), (delta_a, 0, 3*delta_a), coord_system='data')
            proj_1.annotate_line((3*delta_a, 0, delta_a), (-3*delta_a, 0, delta_a), coord_system='data')
            proj_1.annotate_line((3*delta_a, 0, -delta_a), (-3*delta_a, 0, -delta_a), coord_system='data')
            proj_1.annotate_line((-3*delta_a, 0, -delta_a), (-3*delta_a, 0, delta_a), coord_system='data')
            proj_1.annotate_line((3*delta_a, 0, -delta_a), (3*delta_a, 0, delta_a), coord_system='data')
            proj_1.annotate_line((delta_a, 0, -3*delta_a), (-delta_a, 0, -3*delta_a), coord_system='data')
            proj_1.annotate_line((delta_a, 0, 3*delta_a), (-delta_a, 0, 3*delta_a), coord_system='data')
            proj_1.set_cmap(('gas','density'), cmap='magma')
            proj_1.flip_horizontal()
            proj_1.swap_axes()
            colorbar_1 = proj_1.plots['density'].cb
            colorbar_1.vmin = vmin * projection_length
            colorbar_1.vmax = vmax * projection_length
            proj_1.set_axes_unit('mm')
            proj_1.set_ylabel(r'$x \,\mathrm{(mm)}$')
            proj_1.set_xlabel(r'$z \,\mathrm{(mm)}$')
            proj_1.set_font_size(28)
            proj_1.set_zlim('density', vmin*projection_length, vmax*projection_length)
            proj_1.hide_colorbar()
            proj_1.save('y_slice.pdf')
            label = r'Density ($\mathrm{g\,cm^{-3}}$)'
            if 0:
                cfp.plot_colorbar(cmap='magma', vmin=vmin, vmax=vmax, label=label, log=True,
                                    orientation='vertical', panels=3, save='cb_location.pdf')

            noisetype = 'half'
            target_snr = [0.01, 0.1, 1, 3, 5, 10, 100]
            labels = ['SNR = 0.01', 'SNR = 0.1', 'SNR = 1', 'SNR = 3', 'SNR = 5', 'SNR = 10', 'SNR = 100']
            box = 'c_box'
            los = 1
            cg_main,_, _ = get_box(void_size, box, True, los)
            snr_image,_ = noiseimage.add_noise(cg_main, target_snr[-1],noisetype=noisetype)
            vmin = snr_image.min()
            vmax = snr_image.max()
            for val,target in enumerate(target_snr):
                snr_image,_ = noiseimage.add_noise(cg_main, target, noisetype=noisetype)
                cfp.plot_map(snr_image, colorbar=None, vmin=vmin, vmax=vmax)
                plt.axis('off')
                cfp.plot(text=labels[val], normalised_coords=True, x=0.35, y=-0.08, fontsize=18)
                cfp.plot_map(save='./sort_plots/c_box_%d.pdf'%val)
            if 0:
                label = r'Density ($\mathrm{g\,cm^{-2}}$)'
                cfp.plot_colorbar(cmap='magma', vmin=vmin, vmax=vmax, label=label, log=True, orientation='horizontal', panels=2.5, save='cb_snr.pdf')

    timer.report()
