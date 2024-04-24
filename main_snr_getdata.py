
from startlaser import * 
from collections import defaultdict


def get_box(voids, box, projecting, proj_axis=0):
    '''
    Build the box.
    '''
    the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_%s_bubbles/NIF_hdf5_plt_cnt_0375'%(voids))
    the_level = 0
    delta_a = 0.02  #(1.2/3)/2
    edge_mtop = 0.3675 
    edge_mbot = 0.3275

    if box == 'c_box':
        the_ledge = the_ds.arr([-delta_a, edge_mbot, -delta_a],'code_length')
        the_redge = the_ds.arr([delta_a, edge_mtop, delta_a],'code_length')

    if box == 'r_box':
        the_ledge = the_ds.arr([delta_a, edge_mbot, -delta_a],'code_length')
        the_redge = the_ds.arr([delta_a+2*delta_a, edge_mtop, delta_a],'code_length')
    if box == 'f_box':
        the_ledge = the_ds.arr([-delta_a, edge_mbot, delta_a],'code_length')
        the_redge = the_ds.arr([delta_a, edge_mtop, delta_a+2*delta_a],'code_length')

    if box == 'l_box':
        the_ledge = the_ds.arr([-delta_a-2*delta_a, edge_mbot, -delta_a],'code_length')
        the_redge = the_ds.arr([-delta_a, edge_mtop, delta_a],'code_length')
    if box == 'b_box':
        the_ledge = the_ds.arr([-delta_a, edge_mbot, -delta_a-2*delta_a],'code_length')
        the_redge = the_ds.arr([delta_a, edge_mtop, -delta_a],'code_length')

    the_center =0.5*the_ledge + 0.5*the_redge
    mindx = the_ds.index.get_smallest_dx()
    size = ((the_redge - the_ledge)/mindx).astype('int') 

    section = extract_s1.extract_section_1(the_ds,the_level,the_ledge,size)
    cg_first, cg_second, _ = project.Project(section, proj_axis, projecting, field_main, field_two)

    return cg_first, cg_second


def get_spectrumcf(the_box, ncmp=1):
    ''' 
    from cfpack code, plus a few edits/questions:
    1. Q: ncmp = 1  #2 or leave as ncmp=1?; if ncmp > 1, there are more things to add from source.
    2. Edit: norm changed from "forward" to "ortho"; value error below; 1.20.0 version installed for 'forward'
    '''
    data = np.copy(the_box)  #this makes ndim be 1 if it doesn't go into: if (ncmp == 1) and (data.shape[0] > 1) - - -
    if (ncmp == 1) and (data.shape[0] > 1):
        data = np.array([data])
    num = np.array(data[0].shape)
    ndim = len(num)
    kays = -num//2
    ke = np.array([np.min([-kays[d], num[d]//2-1]) for d in range(ndim)])
    
    # 4 nested get_c123_coords defintions
    def get_1d_coords(cmin=0, cmax=1, ndim=10, cell_centred=True):
        if cell_centred:
            d = (cmax-cmin) / float(ndim)
            offset = d/2
        else:
            d = (cmax-cmin) / float(ndim-1)
            offset = 0.0
        return np.linspace(cmin+offset, cmax-offset, num=ndim)
    def get_2d_coords(cmin=[0,0], cmax=[1,1], ndim=[10,10], cell_centred=True):
        cmin = np.array(cmin)
        cmax = np.array(cmax)
        ndim = np.array(ndim)
        if cmin.ndim != 1: cmin = [cmin,cmin]
        if cmax.ndim != 1: cmax = [cmax,cmax]
        if ndim.ndim != 1: ndim = [ndim,ndim]
        c0 = get_1d_coords(cmin=cmin[0], cmax=cmax[0], ndim=ndim[0], cell_centred=cell_centred)
        c1 = get_1d_coords(cmin=cmin[1], cmax=cmax[1], ndim=ndim[1], cell_centred=cell_centred)
        return np.array(np.meshgrid(c0, c1, indexing='ij'))
    def get_3d_coords(cmin=[0,0,0], cmax=[1,1,1], ndim=[10,10,10], cell_centred=True):
        cmin = np.array(cmin)
        cmax = np.array(cmax)
        ndim = np.array(ndim)
        if cmin.ndim != 1: cmin = [cmin,cmin,cmin]
        if cmax.ndim != 1: cmax = [cmax,cmax,cmax]
        if ndim.ndim != 1: ndim = [ndim,ndim,ndim]
        c0 = get_1d_coords(cmin=cmin[0], cmax=cmax[0], ndim=ndim[0], cell_centred=cell_centred)
        c1 = get_1d_coords(cmin=cmin[1], cmax=cmax[1], ndim=ndim[1], cell_centred=cell_centred)
        c2 = get_1d_coords(cmin=cmin[2], cmax=cmax[2], ndim=ndim[2], cell_centred=cell_centred)
        return np.array(np.meshgrid(c0, c1, c2, indexing='ij'))
    def get_coords(cmin, cmax, ndim, cell_centred=True):
        if (type(cmin) != list) and (type(cmin) != np.ndarray): print("need list or nump array inputs", error=True)
        cmin = np.array(cmin)
        cmax = np.array(cmax)
        ndim = np.array(ndim)
        if (cmin.shape != cmax.shape) or (cmax.shape != cmax.shape):
            print("Error: cmin, cmax, ndim, all must have the same shape.", error=True)
        if ndim.shape[0] == 1: return np.array(get_1d_coords(cmin[0], cmax[0], ndim[0], cell_centred))
        if ndim.shape[0] == 2: return np.array(get_2d_coords(cmin, cmax, ndim, cell_centred))
        if ndim.shape[0] == 3: return np.array(get_3d_coords(cmin, cmax, ndim, cell_centred))

    k = get_coords(kays, ke, num, cell_centred=False)
    if ndim == 1: k_abs = np.abs(k)
    if ndim  > 1: k_abs = np.sqrt((k**2).sum(axis=0))
    bins = np.arange(np.max(num)//2) - 0.5
    data_ft = []
    for d in range(ncmp):
        data_ft.append(np.fft.fftn(data[d], norm='forward'))  #ValueError: Invalid norm value forward, should be None or "ortho".
        data_ft[d] = np.fft.fftshift(data_ft[d]) 

    data_ft = np.array(data_ft)

    power_tot = (np.abs(data_ft)**2).sum(axis=0)
    #if ncmp > 1:  # then more things need to happen here 

    def get_binned_stats(data, bin_values, bins=None, statistic='mean', **kwargs):
        from scipy.stats import binned_statistic
        binned_stats = binned_statistic(bin_values.flatten(), data.flatten(), bins=bins, statistic=statistic)
        return binned_stats.statistic, binned_stats.bin_edges
    
    spect_tot, bins = get_binned_stats(power_tot, k_abs, bins=bins)
    bin_centers = bins[:-1]+0.5
    integral_factor = bin_centers**(ndim-1)
    if ndim > 1: integral_factor *= np.pi*2*(ndim-1)
    p_tot = spect_tot*integral_factor
    #ret = {'k': bin_centers, 'P_tot': spect_tot*integral_factor}

    dk = 1/data_ft.size # ADDED

    return bin_centers, p_tot, dk 
    

def noise_iterations(iterations, cg_main, cg_sec, lowpassfilter, target, noise, ksnr=None):
    '''
    Adds noise to the signal if needed. Outputs sigma_rho and rsqrt. 
    '''
    collect_rsq = []
    collect_sigma = []
    collect_sigma3D = []
    for many in range(iterations):  

        # ADD NOISE TO THE IMAGE
        if noise == True:
            cg_wnoise, _ = noiseimage.add_noise(cg_main, target)
        if noise == False:
            cg_wnoise = cg_main
         
        # OBTAIN SIGMA N^2
        sigma_sq, _ = calculate.Calculate(field_main, cg_wnoise, cg_sec) 
        
        ''' will delete soon...
        if lowpassfilter == 'no':
            ks, powspec2D, infxyz = get_spectrumcf(cg_wnoise) 
            #ks, powerspec2D, powerspec3D, infxyz = powerspec.powerspec_lowpass(cg_wnoise)
            #four_amps, powerspec2D, powerspec3D, infxyz, ks = powerspec.powerspec(cg_wnoise, dimension, target)    
            #power_spectrum_map, k = turbulence_driving_analysis.get_power_spectrum(cg_wnoise, zero_pad=False, weights=None, mirror=False, normalise=False, show=False)
        if lowpassfilter == 'yes':
            ks, powspec2D, infxyz = get_spectrumcf(cg_wnoise, the_k=ksnr) 
            #ks, powerspec2D, powerspec3D, infxyz = powerspec.powerspec_lowpass(cg_wnoise, the_k=ksnr, lowpass_filter=lowpassfilter)
            #four_amps, powerspec2D, powerspec3D, infxyz, ks = powerspec.powerspec(cg_wnoise, the_k=ksnr, lowpass_filter=lowpassfilter)
            #power_spectrum_map, k = turbulence_driving_analysis.get_power_spectrum(cg_wnoise, zero_pad=False, weights=None, mirror=False, normalise=False, show=False, lowpass=lowpassfilter, the_k=ksnr)
        ''' 
        # OBTAIN P2D AND P3D           
        ks, powspec2D, infxyz = get_spectrumcf(cg_wnoise) 
        powerspec3D = 2 * ks * powspec2D
        # OBTAIN SIGMA RHO OR ...
        if lowpassfilter == 'no':
            brunts_sigma, sigma_rho2d, brunts_rsq, _ = parsevals.parsevals_th(dimension, powspec2D, powerspec3D, infxyz, ks, sigma_sq) 
        if lowpassfilter == 'yes':
            brunts_sigma, sigma_rho2d, brunts_rsq, _ = parsevals.parsevals_th(dimension, powspec2D, powerspec3D, infxyz, ks, sigma_sq, ksnr) 
        collect_rsq = np.append(collect_rsq, brunts_rsq)
        collect_sigma = np.append(collect_sigma, brunts_sigma.item())  # .item() may not be necessary...

        if noise == False:
            break

    the_rsmed = np.median(collect_rsq)
    the_sigmasmed = np.median(collect_sigma)
    return the_rsmed, the_sigmasmed 


def spectras(cg_signal, target):
    '''
    obtain spectras: signal plus noise, and noise only 
    '''
    for many in range(iterations):  
        cg_wnoise, noiseonly = noiseimage.add_noise(cg_signal, target) 

        ks2D, powspec2D, _ = get_spectrumcf(cg_wnoise) 
        if many == 0:
            collect_spectra = [powspec2D]
        if many > 0:
            collect_spectra = np.append(collect_spectra, [powspec2D], axis=0)  

        ks2D, powspec2D_noise, _ = get_spectrumcf(noiseonly) 
        if many == 0:
            collect_spectranoise = [powspec2D_noise]
        if many > 0:
            collect_spectranoise = np.append(collect_spectranoise, [powspec2D_noise], axis=0) 
         
    return ks2D, collect_spectra, collect_spectranoise



if 0:
    # SETUP!
    # GENERAL BOX SETTINGS 
    field_main = 'density'
    field_two = 'cell_volume'
    polytropic_index = 1 # +/-0.1, for 50mu voids

    projection = True 
    dimension = 1
    noise = True
    iterations = 1000

    # BOX OPTIONS: 15 total per one void size.
    void_size = '12p5mu'
    # 3 different box locations for the x los: left, center, right. and 2 more for the z los: front & back 
    if 1:
        boxes = ['l_box', 'c_box', 'r_box', 'b_box', 'f_box'] 
    if 0: #debug
        boxes = ['c_box'] 
    # 3 different projections
    if 1:  
        projection_axis_list = [0, 1, 2]  
    if 0: #debug
        projection_axis_list = [0, 1]


    # GETTING THE INFO WE NEED TO STORE
    datatype = 'spectra'
    target_snr = [inf, 0.01, 0.1, 1, 3, 5, 10, 100]  

    collect_sigmas = defaultdict(list)
    collect_rsqrtsratio = defaultdict(list)
    collect_rsqrt = defaultdict(list)
    collect_kvals = defaultdict(list)
    collect_rsqrttrue = defaultdict(list)
    collect_rsqrtpure = defaultdict(list)
    rsqrt_full =[]
    collect_spectras = defaultdict(list) 
    collect_noisespectras = defaultdict(list)

    # three projection axes
    for eje, los in enumerate(projection_axis_list):
        # five box locations
        for num, boxing in enumerate(boxes):  
            cg_main, cg_vol = get_box(void_size, boxing, projection, los)
            k_array = np.linspace(1, int(cg_main.shape[0]/2), int(cg_main.shape[0]/2))  #start, stop, num 

            # seven SNR
            for val,target in enumerate(target_snr):
                lowpassfilter = 'no'  

                if datatype != 'spectra':
                    rsqrt_med, sigma_med = noise_iterations(iterations, cg_main, cg_vol, lowpassfilter, target, noise)
                    if target == inf:
                        puresignal = rsqrt_med
                        collect_rsqrtpure[target].append(rsqrt_med)   
                    if target != inf:
                        rsqrt_full = np.append(rsqrt_full, rsqrt_med)
                        collect_sigmas[target].append(sigma_med)  
                        collect_rsqrt[target].append(rsqrt_med)   

                    rsqrt_k = []
                    lowpassfilter = 'yes'  
                    if target != inf:
                        for kval in k_array: 
                            if kval == 1: #for at parsevals only
                                continue
                            rsqrt_med_k, _ = noise_iterations(iterations, cg_main, cg_vol, lowpassfilter, target, noise, kval)
                            rsqrt_k = np.append(rsqrt_k, rsqrt_med_k)
                        index = np.argmin(np.abs(rsqrt_k - puresignal))  #rsqrt_med should be the pure signal value here 
                        kcross = k_array[index] 
                        rsqrt_true, _ = noise_iterations(iterations, cg_main, cg_vol, lowpassfilter, target, noise, kcross)
                        rsqrt_ratio = rsqrt_true/rsqrt_med
                        collect_kvals[target].append(kcross)   
                        collect_rsqrttrue[target].append(rsqrt_true) 
                        collect_rsqrtsratio[target].append(rsqrt_ratio) 
                        print('box snr data ready ', target,boxing,los)

                if datatype == 'spectra':  
                    if target != inf:
                        kvals2D, collectspec, collectnoisespec = spectras(cg_main, target)  
                        if los == projection_axis_list[0] and boxing == boxes[0]:
                            collect_spectras[target] = collectspec
                            collect_noisespectras[target] = collectnoisespec
                        else: 
                            collect_spectras[target] = np.append(collect_spectras[target], collectspec, axis=0) 
                            collect_noisespectras[target] = np.append(collect_noisespectras[target], collectnoisespec, axis=0) 

    if datatype == 'spectra':
        data_spectra = [*collect_spectras.values()]
        data_noisespectra = [*collect_noisespectras.values()]
        # get the mean of all spectra
        the_powspec = np.mean(data_spectra, axis=1)  
        the_noisepowspec = np.mean(data_noisespectra, axis=1)  

    # convert the default dict to a vector array 
    data_sigmas = [*collect_sigmas.values()]
    data_rsqrts = [*collect_rsqrt.values()]
    data_kvals = [*collect_kvals.values()]
    data_rsqrttrue = [*collect_rsqrttrue.values()]
    data_rsqrtratio = [*collect_rsqrtsratio.values()]
    data_rsqrtpure = [*collect_rsqrtpure.values()]

    # and store the information 
    if datatype != 'spectra':
        hfivefile = '/home/552/lv9897/sflaser/h5files/boxplot15_i1kforwardk51ratiofix_%s.h5'%void_size
        Fptr = h5py.File(hfivefile,'w') 
        Fptr['sigmas_%s'%void_size] = data_sigmas
        Fptr['rsqrtfull_%s'%void_size] = data_rsqrts
        Fptr['kvals_%s'%void_size] = data_kvals
        Fptr['rsqrttrue_%s'%void_size] = data_rsqrttrue
        Fptr['ratios_%s'%void_size] = data_rsqrtratio
        Fptr['rsqrtpure_%s'%void_size] = data_rsqrtpure
        Fptr.close()
        print('finished. check h5 file')
    if datatype == 'spectra':
        hfivefile = '/home/552/lv9897/sflaser/h5files/powspec14_i1kSNR_noisesignalwnoise_%s.h5'%void_size
        Fptr = h5py.File(hfivefile,'w')
        Fptr['kvals2D_%s'%void_size] = kvals2D
        Fptr['thepowspec_%s'%void_size] = the_powspec 
        Fptr['thenoisepowspec_%s'%void_size] = the_noisepowspec 
        Fptr.close()
        print('finished. check h5 file')
