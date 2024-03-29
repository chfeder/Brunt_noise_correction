
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


def get_spectrumcf(the_box, ncmp=1, the_k=0):
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
    # added lowpass filter
    if the_k != 0:
        data_ft[d][k_abs > the_k] = 0.0

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

        # OBTAIN P2D AND P3D           
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
        
        # OBTAIN SIGMA RHO OR ...
        powerspec3D = 2 * ks * powspec2D
        brunts_sigma, sigma_rho2d, brunts_rsq, _ = parsevals.parsevals_th(dimension, powspec2D, powerspec3D, infxyz, ks, sigma_sq) 
        collect_rsq = np.append(collect_rsq, brunts_rsq)
        collect_sigma = np.append(collect_sigma, brunts_sigma.item())  # .item() may not be necessary...

        if noise == False:
            break

    the_rsmed = np.median(collect_rsq)
    the_sigmasmed = np.median(collect_sigma)
    return the_rsmed, the_sigmasmed 


def spectras(cg_signal, target, justnoise=False):
    '''
    obtain spectras: signal plus noise, and noise only 
    '''
    for many in range(iterations):  
        cg_wnoise, noiseonly = noiseimage.add_noise(cg_main, target) #cg_signal CHECK
        if justnoise == False:
            ks2D, powspec2D, _ = get_spectrumcf(cg_wnoise) 
            if many == 0:
                collect_spectra = [powspec2D]
            if many > 0:
                collect_spectra = np.append(collect_spectra, [powspec2D], axis=0) 
            collect_spectranoise = 0

        if justnoise == True:
            print('target in spectras def ', target)
            ks2D, powspec2D_noise, _ = get_spectrumcf(noiseonly) 
            if many == 0:
                collect_spectranoise = [powspec2D_noise]
            if many > 0:
                collect_spectranoise = np.append(collect_spectranoise, [powspec2D_noise], axis=0) 
            collect_spectra = 0
         
    return ks2D, collect_spectra, collect_spectranoise




# SETUP!
# GENERAL BOX SETTINGS 
field_main = 'density'
field_two = 'cell_volume'
polytropic_index = 1 # +/-0.1, for 50mu voids

projection = True 
dimension = 1
noise = True
iterations = 2

# for pure signal, simply do target_snr[-1], for all others target_snr[:-1]
target_snr = [0.01, 0.1, 1, 3, 5, 10, 100, inf]  #edit to see effect on 0.001 and 0.0001

# BOX OPTIONS: 15 total per one void size.
void_size = '50mu'
# 3 different box locations for the x los: left, center, right. and 2 more for the z los: front & back 
if 0:
    boxes = ['l_box', 'c_box', 'r_box', 'b_box', 'f_box'] 
if 1: #debug
    boxes = ['c_box'] 
# 3 different projections
if 0:  
    projection_axis_list = [0, 1, 2]  
if 1: #debug
    projection_axis_list = [0, 1]




# GETTING THE INFO WE NEED TO STORE
datatype = 'allother'
if 1:
    collect_sigmas = defaultdict(list)
    collect_rsqrtsratio = defaultdict(list)
    collect_rsqrt = defaultdict(list)
    collect_kvals = defaultdict(list)
    collect_rsqrttrue = defaultdict(list)
    collect_rsqrtpure = defaultdict(list)
    rsqrt_full =[]
    collect_spectras = defaultdict(list) 
    collect_noisespectra = defaultdict(list)

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
                    if target != inf:
                        rsqrt_full = np.append(rsqrt_full, rsqrt_med)
                        collect_sigmas[target].append(sigma_med)  
                        collect_rsqrt[target].append(rsqrt_med)   
                    if target == inf:
                        collect_rsqrtpure[target].append(rsqrt_med)   

                if datatype == 'spectra':
                    kvals2D, collectspec, _ = spectras(cg_main, target)  
                    if target != inf:
                        if los == projection_axis_list[0] and boxing == boxes[0]:
                            collect_spectras[target] = collectspec
                        else: 
                            collect_spectras[target] = np.append(collect_spectras[target], collectspec, axis=0) 

                        if los == projection_axis_list[-1] and boxing == boxes[-1]:
                            justnoise = True
                            print('target in execution ', target)
                            _, _, collectnoisespec = spectras(cg_main, target, justnoise)                  
                            collect_noisespectra[target] = collectnoisespec

            # could have done all at once, but found it less tedious to have rsqrt full first  
            if datatype != 'spectra':
                for val,target in enumerate(target_snr):
                    rsqrt_k = []
                    lowpassfilter = 'yes'  
                    if target != inf:
                        for kval in k_array: 
                            rsqrt_med_k, _ = noise_iterations(iterations, cg_main, cg_vol, lowpassfilter, target, noise, kval)
                            rsqrt_k = np.append(rsqrt_k, rsqrt_med_k)
                        index = np.argmin(np.abs(rsqrt_k - rsqrt_med))  #rsqrt_med should be the pure signal value here 
                        kcross = k_array[index] 
                        rsqrt_true, _ = noise_iterations(iterations, cg_main, cg_vol, lowpassfilter, target, noise, kcross)
                        rsqrt_ratio = rsqrt_true/rsqrt_full[val]
                        collect_kvals[target].append(kcross)   
                        collect_rsqrttrue[target].append(rsqrt_true) 
                        collect_rsqrtsratio[target].append(rsqrt_ratio) 
                        print('box snr data ready ', target,boxing,los)

    if datatype == 'spectra':
        data_spectra = [*collect_spectras.values()]
        data_noisespectra = [*collect_noisespectra.values()]
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
        hfivefile = '/home/552/lv9897/sflaser/h5files/boxplot15_i2forwardk51_newpipe_%s.h5'%void_size
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
        hfivefile = '/home/552/lv9897/sflaser/h5files/powspec14_i2SNRforward_newpipe_%s.h5'%void_size
        Fptr = h5py.File(hfivefile,'w')
        Fptr['kvals2D_%s'%void_size] = kvals2D
        Fptr['thepowspec_%s'%void_size] = the_powspec 
        Fptr['thenoisepowspec_%s'%void_size] = the_noisepowspec 
        Fptr.close()
        print('finished. check h5 file')




# GETTING THE INFO FROM STORAGE
if 0:
    # INITIALIZE THE FIGURE
    fig,ax = plt.subplots(1,1)
    voids = ['50mu']#'12p5mu, '50mu', '100mu']
    colors = ['g'] #'orange','g', 'b']
    for draw, void_size in enumerate(voids):
        whichfig = 'sigmas_%s'%voids[0]
        whichfigtwo = None 
        
        #hfivefile = '/home/552/lv9897/sflaser/h5files/powspec14_i1k_newpipe_%s.h5'%void_size
        hfivefile = '/home/552/lv9897/sflaser/h5files/boxplot15_i2forwardk51simplify_newpipe_%s.h5'%void_size
        Fptr = h5py.File(hfivefile,'r')
        data = list(Fptr['%s'%whichfig])  #returns a list 
        #stop()
        if whichfigtwo:
            datatwo = list(Fptr['%s'%whichfigtwo])    


        # POWERSPECTRA 
        if whichfig == 'kvals2D':
            whichfigtwo = 'thepowspec'
            whichfigthree = 'thenoisepowspec'
            xdata = Fptr['%s'%whichfig][()]  #returns a numpy array
            if whichfigtwo:
                ydata = Fptr['%s'%whichfigtwo][()]  
            if whichfigthree:
                nydata = Fptr['%s'%whichfigthree][()]  

            for i in range(len(target_snr[:-1])):
                ax.plot(xdata, ydata[i], color='blue', alpha=0.5, label='powspec') 
                ax.plot(xdata, nydata[i], color='grey', label='powspec') 
            ax.set(xscale='log', yscale= 'log', xlabel='k values', ylabel=r'$P_{2D}$')#, xlim=(5,51)) #, ylim=(10e-1,10e11))#, xlim=(1e0,1e2))  #xscale= 'log',
            #ax.axvline(x=51)
            plt.savefig('sort_plots/powspec14_i2SNRforward_newpipe_%s'%(void_size))
            print('saved spectra14!')
            plt.clf()
            Fptr.close()         

        # BOXPLOTS WITH SCATTERS
        if whichfig != 'kvals2D':
            xdata = np.array(target_snr[:-1])
            xdatalog = np.log10(xdata)  
            # funny order to make the boxplots with plot overlay to work
            plt.xscale('log')
            plt.xlim(np.array(xdata).min()/1.20, np.array(xdata).max()*1.20)

            # the boxplot!  WARNING: this sets only for 'sigmas' option 
            the_data = data[-1]/data
            data_list = the_data.tolist()
            boxplot = ax.boxplot(data, showfliers=True, showmeans=False, labels=xdata, positions=xdata)#, boxprops=dict(color="none"), capprops=dict(color="none"), medianprops=dict(color="none"))  

            # ...and hacking the box plot 
            cap_s = 1.08
            box_s = 1.16
            for i in range(len(data)):
                cap = [xdata[i]/cap_s, xdata[i]*cap_s]
                boxplot['caps'][2*i+0].set_data([cap, boxplot['caps'][2*i+0].get_data()[1]])  #lower cap
                boxplot['caps'][2*i+1].set_data([cap, boxplot['caps'][2*i+1].get_data()[1]])  #upper cap

                med = [xdata[i]/box_s, xdata[i]*box_s]
                boxplot['medians'][i].set_data([med, boxplot['medians'][i].get_data()[1]]) 
                
                box = [xdata[i]/box_s, xdata[i]*box_s, xdata[i]*box_s, xdata[i]/box_s, xdata[i]/box_s]
                boxplot['boxes'][i].set_data([box, boxplot['boxes'][i].get_data()[1]])

            medians = [median.get_ydata()[0] for median in boxplot['medians']]
            #ydatalog = np.log10(medians)
            #stop()
            plt.scatter(xdata, medians, color=colors[draw], marker='o', label='Median', alpha=0.5)  #reduce the marker size! 

            if whichfigtwo:
                medians_two=[]
                # FOR 15 boxes, all target
                if 1:
                    for k in range(7):
                        medians_two = np.append(medians_two, np.median(datatwo[k]))
                    plt.scatter(xdata, medians_two, color='orange', marker='o', alpha=0.5, label='other median')  #reduce the marker size! 
                # FOR 15 boxes, one target - signal only
                if 0:
                    medians_two = np.array([np.median(datatwo)]*len(xdata))
                    plt.scatter(xdata, medians_two, color='orange', marker='o', alpha=0.5, label='other median')  #reduce the marker size! 


            # POSSIBLE FIT FUNCTIONS
            fitfunct = 'linear'
            def linfunct(x, a, b):  #I may need to take the log of this!!
                y = a * x + b
                return y
            def tanfunct(x, a, b, c, d):
                y = a * np.tanh((c*x)+d) + b
                return y

            def lininvfunct(x, a, b=1):
                y = a * x**(-1) + b
                return y
            def tanfunct_simple(x, first, second):  #used to be x, c, d, e; try this first when data is ready
                y = first * np.tanh(second*np.log10(x)) + (1-first)  #c should be negative, d positive; c+e=1
                return y

            def sigmarho_fit(x,a):#, c,d,e):
                if void_size=='12p5mu':
                    first =1.09314738 
                    second = -1.8465386429
                    third = 2.039733348
                if void_size=='50mu':
                    first =0.990295945
                    second =-2.0681085568
                    third =1.9667865114
                if void_size=='100mu':
                    first = 1.0847199980 
                    second = -2.1120163828 
                    third = 2.14165109959
                y = 1/ (lininvfunct(x, a, b=1.0) * tanfunct_simple(x, first,second,third)) 
                return y 

            def sigmfunct(x, a, b):
                y = 1/(1 + a*np.exp(-x+b)) 
                return y

            # GET the median, the 16th ('-' correction), and the 84th ('+' correction) of fit params a, b, c, d
            # this then allows us to report param +/- correction 
            if 0:
                param_a = []
                param_b = []
                param_c = []
                param_d = []
                param_e = []
                for j in range(1):  #avoid hard-coding 15
                    oneboxproj = []
                    for k in range(len(target_snr[:-1])):
                        prexdata = data[k][j]
                        prexdatalog = np.log10(prexdata)
                        if 0:
                            oneboxproj = np.append(oneboxproj, prexdatalog) 
                        if 1:
                            oneboxproj = np.append(oneboxproj, prexdata) 
                            #oneboxproj = np.append(oneboxproj, data[-1][j]/prexdata)  #sigma_rho(SNR->inf)/sigma_rho(SNR) 
                    params = {"first":[-2,-1,-0.1], "second":[0.5,2.0,4]} #h, "d":[-3,-2.5,-2], "e":[1.75,2,2.25]}
                    #params = {"a":[-1.3,0.0,1.35], "b":[-2.15,0,2.15], "c":[1.60,2.1,2.60]} 
                    #rets = cfp.fit(sigmarho_fit, xdata, oneboxproj, params=params) 
                    rets = cfp.fit(tanfunct_simple, xdata, oneboxproj, params=params) 
                        
                    param_a = np.append(param_a, rets.popt[0])
                    if 1:
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
                if 1:
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

                if 0: 
                    med_b = p16_b = p84_b = 0
                if 1: 
                    med_c = p16_c = p84_c = 0
                if 1:
                    med_d = p16_d = p84_d = 0
                if 1:
                    med_e = p16_e = p84_e = 0

                param_meds = [med_a, med_b]#, med_c, med_d, med_e]
                print("param meds ",*param_meds)
                stop()
                param_16s = [p16_a, p16_b]#, p16_c, p16_d, p16_e]  #report these as - for param
                param_84s = [p84_a, p84_b]#, p84_c, p84_d, p84_e]  #report these as + for param

            the_xrecipe = np.linspace(xdata.min(),xdata.max(),num=500)  #2000  
            #the_xten = 10**the_xrecipe
            the_xten = the_xrecipe

            # fitting for the params of the boxplot
            # to fit the median explicitly with the params..
            if 0:
                the_yrecipemed = tanfunct_simple(the_xrecipe, *param_meds)  
                the_ytenmed = the_yrecipemed 

                the_yrecipelow = tanfunct_simple(the_xrecipe, *param_16s)
                the_ytenlow = the_yrecipelow

                the_yrecipeupp = tanfunct_simple(the_xrecipe, *param_84s)
                the_ytenupp = the_yrecipeupp
            
            if 0:
                # ..OR with the median of the data
                ydata_med = np.percentile(data, 50, axis=1)  #axis=0 for the 15 boxes, axis=1 for the SNRs 
                ydatalog_med = np.log10(ydata_med)
                if 0:
                    rets_med = cfp.fit(lininvfunct, xdata, ydatalog_med)  
                if 1:
                    rets_med = cfp.fit(sigmarho_fit, xdata, ydata_med[-1]/ydata_med)  
                    #rets_med = cfp.fit(tanfunct_simple, xdata, ydata_med)  
                the_yrecipemed = sigmarho_fit(the_xrecipe, *rets_med.popt)
                #the_yrecipemed = tanfunct_simple(the_xrecipe, *rets_med.popt)
                if 0:
                    the_ytenmed = 10**the_yrecipemed
                if 1:
                    the_ytenmed = the_yrecipemed

                if 1:
                    # to fit the lower bound, 16th
                    ydata_low = np.percentile(data, 16, axis=1) 
                    rets_low = cfp.fit(sigmarho_fit, xdata, ydata_low[-1]/ydata_low)
                    the_yrecipelow = sigmarho_fit(the_xrecipe, *rets_low.popt)
                    the_ytenlow = the_yrecipelow

                    # to fit the upper bound, 84th
                    ydata_upp = np.percentile(data, 84, axis=1) 
                    rets_upp = cfp.fit(sigmarho_fit, xdata, ydata_upp[-1]/ydata_upp)
                    the_yrecipeupp = sigmarho_fit(the_xrecipe, *rets_upp.popt)
                    the_ytenupp = the_yrecipeupp 

            # ...and plot!
            if 0:
                plt.plot(the_xten, the_ytenmed, color=colors[draw], linestyle='dashed')
            if 0:
                plt.plot(the_xten, the_ytenlow, color=colors[draw], linestyle='dotted', alpha=0.7, label='lower') 
                plt.plot(the_xten, the_ytenupp, color=colors[draw], linestyle='dotted', alpha=0.7, label='upper')  
                plt.fill_between(the_xten, the_ytenlow, the_ytenupp, color=colors[draw], alpha=0.2)

            if 0:  #pass the params needed; not a fit
                first =0.990295945  #name these with more sensible names; saturation, exponent, etc.
                second =-2.0681085568
                third =1.9667865114
                fixed = 1/(third+first)
                print("%s"%void_size)
                plt.plot(the_xten, sigmarho_fit(the_xten, fixed, first, second, third), color=colors[draw],linestyle='dashed', alpha=0.4) 

    # SAVE    
    plt.yscale('log')
    plt.xlabel('SNR')

    # all the labels, change full to noise
    #plt.ylabel(r'$\sigma_{\rho/\rho_o,inf} / \sigma_{\rho/\rho_o,snr}$') 
    plt.ylabel(r'$\sigma_{\rho/\rho_o ,Br}$') 
    #plt.ylabel(r'$R^{1/2}_{noisy}$')
    #plt.ylabel('k values')
    #plt.ylabel(r'$R^{1/2}_{denoised}$, $R^{1/2}_{pure signal}$')
    #plt.ylabel(r'$R^{1/2}_{noise}$, $R^{1/2}_{denoised}$')
    #plt.ylabel(r'$R^{1/2}_{denoised}/R^{1/2}_{noisy}$')

    # really saving now
    plt.savefig('sort_plots/boxplot15_%s_i2newpipe51'%whichfig)
    #plt.savefig('sort_plots/boxplot15_%s_ratio_abfreenorhalf_noboxes_newpipe'%(whichfig))
    print('%s_saved!'%whichfig)
    plt.clf()
    Fptr.close()




## sigma rho vs snr multiplied, all free:
# 50um:  a =  2.1579298, b =  0.39977052868, c =  0.0174629205, d =  -44.6403574, e =  0.4143853683
## sigma rho vs snr divide: 
# 50um: a =  19.8668843, b =  17.3482837809, c =  -18.8323664, d =  -1.4965170289, e =  40.1989922

## y = a * np.tanh(b*log(x)) + c  
# 12p5um: a= 1.09314738, b = -1.8465386429, c = 2.039733348
# 50um:  a =  0.990295945, b =  -2.0681085568, c =  1.9667865114371 
# 100um: a =  1.0847199980, b= -2.1120163828, c=2.14165109959

# Sigma_B vs SNR linear fits:
# y = a * x**(-1) + b
# 12.5um: a=same-sh, b=0.1233595; bfix=0.2, a= 0.8831989 bfix=0.3, a=0.882093696; bfix=0.25, a=0.8826463229601 
# 50.0um: a=0.93244, b = 0.154046, bfixed=0.3: a=0.9304810
# 100.um: a=same-sh, b=0.2279368; bfixed=0.4, a=1.7156249314, bfixed=0.5(hits 100 exactly), a=1.714519679

# PLAN 
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
