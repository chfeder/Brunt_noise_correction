
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





# GETTING THE INFO FROM STORAGE: NEXT, separate things into definitions! 
if 1:
    # POSSIBLE FIT FUNCTIONS
    def linfunct(x, a, b):  
        y = a * x + b
        return y

    def lininvfunct(x, a, b=1):
        y = a * x**(-1) + b
        return y

    def tanfunct_simple(x, first, second):  #used to be x, c, d, e; try this first when data is ready
        y = first * np.tanh(second*np.log10(x)) + (1-first)  #c should be negative, d positive; c+e=1
        return y

    def sigmarho_fit(x,a):  
        y = 1/ (lininvfunct(x, a, b=1.0) * tanfunct_simple(x, first=first, second=second)) 
        return y 


    def fitthebox(whichfigure, eldata, secondparam):  
        '''
        GET the median, the 16th ('-' correction), and the 84th ('+' correction) of fit params a, b, c, d, etc
        this then allows us to report param +/- correction 
        '''
        param_a = []
        param_b = []
        param_c = []
        param_d = []
        param_e = []
        for j in range(15):  #avoid hard-coding numbers, but there are 15 boxes so...
            oneboxproj = []
            for k in range(len(xdata)):
                prexdata = eldata[k][j]
                if whichfigure == 'sigmas':
                    oneboxproj = np.append(oneboxproj, prexdata)  #sigma_rho(SNR->inf)/sigma_rho(SNR) 
                    #oneboxproj = np.append(oneboxproj, eldata[-1][j]/prexdata)  #sigma_rho(SNR->inf)/sigma_rho(SNR) 
                else: 
                    oneboxproj = np.append(oneboxproj, prexdata) 
            
            if whichfigure == 'sigmas':
                params=None
                rets = cfp.fit(sigmarho_fit, xdata, oneboxproj, params=params) 
            else: 
                params = {"first":[-2,-1,-0.1], "second":[0.5,2.0,4]} 
                rets = cfp.fit(tanfunct_simple, xdata, oneboxproj, params=params) 
                
            param_a = np.append(param_a, rets.popt[0])
            if secondparam == 1:
                param_b = np.append(param_b, rets.popt[1])
            if 0:
                param_c = np.append(param_c, rets.popt[2])
            if 0:
                param_d = np.append(param_d, rets.popt[3])
            if 0:
                param_e = np.append(param_e, rets.popt[4])

        med_a = np.percentile(param_a, 50) #same as np.median()
        p16_a = np.percentile(param_a, 16)
        p84_a = np.percentile(param_a, 84)
        if secondparam == 1:
            med_b = np.percentile(param_b, 50)
            p16_b = np.percentile(param_b, 16)
            p84_b = np.percentile(param_b, 84)
        if 0:
            med_c = np.percentile(param_c, 50)
            p16_c = np.percentile(param_c, 16)
            p84_c = np.percentile(param_c, 84)
        if 0:
            med_d = np.percentile(param_d, 50)
            p16_d = np.percentile(param_d, 16)
            p84_d = np.percentile(param_d, 84)
        if 0:
            med_e = np.percentile(param_e, 50)
            p16_e = np.percentile(param_e, 16)
            p84_e = np.percentile(param_e, 84)

        if secondparam == 0:
            med_b = p16_b = p84_b = 0
        if 1: 
            med_c = p16_c = p84_c = 0
        if 1:
            med_d = p16_d = p84_d = 0
        if 1:
            med_e = p16_e = p84_e = 0

        if secondparam == 0:
            param_meds = [med_a]
            param_16s = [p16_a] 
            param_84s = [p84_a] 
        else:
            param_meds = [med_a, med_b]
            param_16s = [p16_a, p16_b]#, p16_c, p16_d, p16_e]  #report these as - for param
            param_84s = [p84_a, p84_b]#, p84_c, p84_d, p84_e]  #report these as + for param
            #stop()
        
        return param_meds, param_16s, param_84s


    def plotthefitting(thefigures, themedians, sixteenths, eightyfours):
        if thefigures == 'sigmas':
            the_yrecipemed = sigmarho_fit(the_xrecipe, *themedians)  
            the_yrecipelow = sigmarho_fit(the_xrecipe, *sixteenths)
            the_yrecipeupp = sigmarho_fit(the_xrecipe, *eightyfours)
        if thefigures == 'ratios':
            the_yrecipemed = tanfunct_simple(the_xrecipe, *themedians)  
            the_yrecipelow = tanfunct_simple(the_xrecipe, *sixteenths)
            the_yrecipeupp = tanfunct_simple(the_xrecipe, *eightyfours)
        return the_yrecipemed, the_yrecipelow, the_yrecipeupp


    def boxplotscatter(thefigure, ladata, thesecond=0, secondfig=None, datasecond=None):
        fig,ax = plt.subplots(1,1)

        # funny order to make the boxplots with a plot overlayed to work
        plt.xscale('log')
        plt.xlim(np.array(xdata).min()/1.20, np.array(xdata).max()*1.20)

        # the boxplot!  WARNING: this sets only for 'sigmas' option 
        if thefigure == 'sigmas':
            the_data = ladata[-1]/ladata
            data_lista = the_data.tolist()
        else:
            data_lista = ladata
        boxplot = ax.boxplot(data_lista, showfliers=True, showmeans=False, labels=xdata, positions=xdata) #\
                            #, boxprops=dict(color="none"), capprops=dict(color="none"), medianprops=dict(color="none"))  

        # ...and hacking the box plot 
        cap_s = 1.08
        box_s = 1.16
        for i in range(len(ladata)):
            cap = [xdata[i]/cap_s, xdata[i]*cap_s]
            boxplot['caps'][2*i+0].set_data([cap, boxplot['caps'][2*i+0].get_data()[1]])  #lower cap
            boxplot['caps'][2*i+1].set_data([cap, boxplot['caps'][2*i+1].get_data()[1]])  #upper cap

            med = [xdata[i]/box_s, xdata[i]*box_s]
            boxplot['medians'][i].set_data([med, boxplot['medians'][i].get_data()[1]]) 
            
            box = [xdata[i]/box_s, xdata[i]*box_s, xdata[i]*box_s, xdata[i]/box_s, xdata[i]/box_s]
            boxplot['boxes'][i].set_data([box, boxplot['boxes'][i].get_data()[1]])

        medians = [median.get_ydata()[0] for median in boxplot['medians']]
        plt.scatter(xdata, medians, color=colors[draw], marker='o', label='Median', alpha=0.5)  #reduce the marker size! 

        # for overlaying purposes, such as rhalf noise, rhalf denoised; may end up not needing/otherwise revise/edit
        if secondfig: 
            medians_two=[]
            # FOR 15 boxes, all target
            if 0:
                for k in range(7):
                    medians_two = np.append(medians_two, np.median(datasecond[k]))
                plt.scatter(xdata, medians_two, color='orange', marker='o', alpha=0.5, label='other median')  #reduce the marker size! 
            # FOR 15 boxes, one target - signal only
            if 0:
                medians_two = np.array([np.median(datasecond)]*len(xdata))
                plt.scatter(xdata, medians_two, color='orange', marker='o', alpha=0.5, label='other median')  #reduce the marker size! 


        # call the fittings procedure definition
        meds, lows, upps = fitthebox(thefigure, data_lista, thesecond)
        # call the fitting definition, which will provide the the_yrecipes below
        if thefigure =='ratios':
            ymed, ylow, yupp = plotthefitting(thefigure, meds, lows, upps)

            plt.plot(the_xrecipe, ymed, color=colors[draw], linestyle='dashed')
            plt.plot(the_xrecipe, ylow, color=colors[draw], linestyle='dotted', alpha=0.7, label='lower') 
            plt.plot(the_xrecipe, yupp, color=colors[draw], linestyle='dotted', alpha=0.7, label='upper')  
            plt.fill_between(the_xrecipe, ylow, yupp, color=colors[draw], alpha=0.2)

        # SAVE    
        # ALSO AUTOMATE SOOOOOON.
        #plt.ylabel(r'$\sigma_{\rho/\rho_o ,Br}$') 
        #plt.ylabel(r'$R^{1/2}_{noisy}$')
        #plt.ylabel(r'$R^{1/2}_{denoised}$, $R^{1/2}_{pure signal}$')
        #plt.ylabel(r'$R^{1/2}_{noise}$, $R^{1/2}_{denoised}$')
        if thefigure == 'kvals':
            plt.yscale('log')
            plt.ylabel('k noise')
        if thefigure == 'sigmas':
            plt.yscale('log')
            plt.ylabel(r'$\sigma_{\rho/\rho_o,inf} / \sigma_{\rho/\rho_o,snr}$') 
        if thefigure == 'ratios':
            plt.ylabel(r'$R^{1/2}_{denoised}/R^{1/2}_{noisy}$')
        plt.xlabel('SNR')

        # really saving now
        plt.savefig('sort_plots/boxplot15_%s_i1kfix_%s'%(thefigure,void_size))
        print('%s_saved!'%thefigure)
        plt.clf()
        Fptr.close() 
        return meds



    # LOOP FOR THREE VOIDS!
    target_snr = [0.01, 0.1, 1, 3, 5, 10, 100, inf]  
    xdata = np.array(target_snr[:-1])
    the_xrecipe = np.linspace(xdata.min(), xdata.max(), num=2500)  
    voids = ['50mu', '100mu', '12p5mu']
    colors = ['g', 'b', 'orange']
    for draw, void_size in enumerate(voids):
        print('void size is ',void_size)
        whichfig = 'kvals'
        whichfigtwo = None 
        mainfigure = None #'sigmas'


        # ALL OTHER DATA: BOXPLOTS + SCATTERS + FITS
        if whichfig != 'kvals2D':
            hfivefile = '/home/552/lv9897/sflaser/h5files/boxplot15_i1kforwardk51ratiofix_%s.h5'%void_size
            Fptr = h5py.File(hfivefile,'r')
            # these are parameters for the tanfunct_simple fit, initialized to zero for multiple plot purposes 
            first = 0
            second = 0

            data = list(Fptr['%s_%s'%(whichfig,void_size)])  #returns a list 
            if whichfigtwo:
                data_two = list(Fptr['%s_%s'%(mainfigure,void_size)])  #returns a list 

            if whichfig == 'kvals':
                for val, kvalues in enumerate(data):
                    data[val] = data[val]/51 
            boxplotscatter(whichfig, data, thesecond=1)  #see if it still works.
            if mainfigure == 'sigmas':
                first = listofmeds[0]
                second = listofmeds[1]
                boxplotscatter(mainfigure, data_two, thesecond=0)


        # POWERSPECTRA: #make this a small definition too?? 
        if whichfig == 'kvals2D':
            hfivefile = '/home/552/lv9897/sflaser/h5files/powspec14_i1kSNR_noisesignalwnoise_%s.h5'%void_size
            Fptr = h5py.File(hfivefile,'r')
            fig,ax = plt.subplots(1,1)
            
            whichfigtwo = 'thepowspec'
            whichfigthree = 'thenoisepowspec'
            ps_xdata = Fptr['%s_%s'%(whichfig,void_size)][()]  #returns a numpy array
            if whichfigtwo:
                ydata = Fptr['%s_%s'%(whichfigtwo,void_size)][()]  
            if whichfigthree:
                nydata = Fptr['%s_%s'%(whichfigthree,void_size)][()]  

            for i in range(len(target_snr[:-1])):
                ax.plot(ps_xdata, ydata[i], color='blue', alpha=0.5, label='powspec') 
                ax.plot(ps_xdata, nydata[i], color='grey', label='powspec') 
            ax.set(yscale= 'log', xlabel='k values', ylabel=r'$P_{2D}$')#, xlim=(5,51)) #, ylim=(10e-1,10e11))#, xlim=(1e0,1e2))  #xscale= 'log',
            #ax.axvline(x=51)
            plt.savefig('sort_plots/powspec14_i1kSNRforward_noisesignalwnoise_%s'%(void_size))
            print('saved spectra14!')
            plt.clf()
            Fptr.close()         



# IDAL NEXT: make this a 'signal to noise ratio study' CLASS and call it from a truer main.py file or from the bottom of this class? 


# PLAN : delete from here soon.
'''
i. slice of positions of the cubes 
ii. plot 2d power spectrums of noise only and signal + noise

# for any box(position, proj, void) for each SNR 
a. get R^{1/2}_full values, where final is R^{1/2}_pure - call SNR function(k values turned off)

b. find the k values: AUTOMATE [DONE] 
    1. based on ii, for a possible range of k values for a specific SNR target, boxplot R^{1/2} vs that k range; iterations = 1000
    2. pick the k value that matches R^{1/2} with R^{1/2}_inf most closely BY HAND; store 
# pseudocode on automating
index = np.argmin(np.abs( r^1/2array - r1/2value ))
kcross = k_crossarray[index]
# these k values to recover rsqrt_true, these were picked by hand, could now check matching 
#ksnr = [6, 8, 18, 29, 35, 44, 52]  #right, left:[5, 7, 16, 25, 30, 38, 52], center: ksnr = [8, 8, 9, 19, 30, 36, 44, 52]

c. get R^{1/2}_true values - call SNR function(k values turned on)

# figures for any box(position, proj, void) or boxes combined 'superplot'
# plot or boxplot ...
1. R^{1/2}_full vs SNR
2. R^{1/2}_true vs SNR with R^{1/2}_pure
3a. R^{1/2}_true/R^{1/2}_full vs SNR
3b. Fit 3a. - Make a table with all the parameters
4?. k values vs SNR
'''

'''
if 0:  #working on the cfp tool
y_err_up = the_ytenupp - the_ytenmed 
y_err_low = the_ytenmed - the_ytenlow 
y_err_arr = [y_err_low,y_err_up] 
cfp.plot(x=the_xten, y=the_ytenmed, yerr=y_err_arr, label='params median', shaded_err=True, save='testingthis.png', ylog=True)
'''
