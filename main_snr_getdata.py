
from startlaser import * 
from collections import defaultdict


def get_box(voids, box, projecting, proj_axis=0, imaging='no'):
    '''
    Build the box. Get the cover grid. Project (or not)
    '''
    the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_512_%s_bubbles/NIF_hdf5_plt_cnt_0375'%(voids))
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
    cg_first, cg_first_norm, cg_second = project.Project(section, proj_axis, projecting, field_main, field_two)

    return cg_first, cg_first_norm, cg_second


def get_spectrumcf(the_box, ncmp=1, sum=False, return_ft_data=True, verbose=False): 

    # check statistic type
    if sum:
        statistic = 'sum'
        normalise = False
    else:
        statistic = 'mean'
        normalise = True

    # flatten to strip any possible extra dimensions
    data = np.copy(the_box)  
    if (ncmp == 1) and (data.shape[0] > 1):
        data = np.array([data])  # add an extra (fake) index, so we can index as if there were components
    num = np.array(data[0].shape)  # number of points in data
    ndim = len(num)  # number of dimensions

    # print info about real-space field
    mean_sq = 0.0
    std_sq = 0.0
    for d in range(ncmp):
        mean_sq += data[d].mean()**2
        std_sq += np.std(data[d])**2
    if verbose: print("data mean squared = "+cfp.eform(mean_sq)+", data std squared = "+cfp.eform(std_sq))

    # construct wave numbers and bins
    ks = -(num//2)  # k start
    ke = np.array([num[d]//2+(num[d]%2-1) for d in range(ndim)]) # k end
    k = cfp.get_coords(ks, ke, num, cell_centred=False) # get k vector with k=0 in the center
    if ndim == 1: k_abs = np.abs(k) # length of k vector
    if ndim  > 1: k_abs = np.sqrt((k**2).sum(axis=0)) # length of k vector
    bins = np.arange(np.max(ke)+2) - 0.5 # k bins for 1D spectrum

    # do Fourier transformation
    data_ft = []
    for d in range(ncmp):
        data_ft.append(np.fft.fftn(data[d].T, norm='forward'))  # FFT (note that we transpose before FFT, to get data[x,[y,[z]]] -> data[[[z],y],x])
        data_ft[d] = np.fft.fftshift(data_ft[d].T)  # shift k=0 to center and transpose back, so we have data_ft[kx,[ky,[kz]]]
    data_ft = np.array(data_ft)
    
    # get total power
    power_tot = (np.abs(data_ft)**2).sum(axis=0)
    P0 = power_tot[k_abs==0][0] # extract k=0 value of power spectrum
    norm = power_tot.sum()-P0 # integral over all k!=0
    if verbose: print("           P(k=0) = "+cfp.eform(P0)+",      total power = "+cfp.eform(norm))

    # bin in k shells
    def get_binned_stats(data, bin_values, bins=None, statistic='mean', **kwargs):
        from scipy.stats import binned_statistic
        binned_stats = binned_statistic(bin_values.flatten(), data.flatten(), bins=bins, statistic=statistic)
        return binned_stats.statistic, binned_stats.bin_edges
    spect_tot, bins = get_binned_stats(power_tot, k_abs, bins=bins, statistic=statistic) 
    bin_centers = bins[:-1]+0.5
    if statistic == "mean":
        integral_factor = bin_centers**(ndim-1)
        if ndim > 1: integral_factor *= np.pi*2*(ndim-1)

    # integral factor and normalisation if needed
    if statistic == "sum": integral_factor = np.ones(spect_tot.shape)
    spect_tot[bin_centers!=0] *= integral_factor[bin_centers!=0]  # k=0 stays
    if normalise: spect_tot[bin_centers!=0] = spect_tot[bin_centers!=0] / spect_tot[bin_centers!=0].sum() * norm # normalise

    
    # construct return dict - if needed 
    ret = {'k': bin_centers, 'P_tot': spect_tot}
    if return_ft_data: ret['Power_tot'] = power_tot # if the user wants to have the Fourier-transformed data (power) returned
    #return ret # return dict with power spectrum data    

    return bin_centers, spect_tot, power_tot 
    

def noise_iterations(iterations, cg_main_sum, cg_sec, lowpassfilter, target, ksnr=None, correct=None, typeofnoise=None):
    collect_rsq = []
    collect_sigma = []
    collect_sigmaN = []

    for many in range(iterations):  
        # ADD NOISE TO THE IMAGE
        cg_wnoise_sum, noiseonly = noiseimage.add_noise(cg_main_sum, target, noisetype=typeofnoise)
        cg_wnoise, _ = noiseimage.add_noise(correct, target, noisetype=typeofnoise)
     
        # OBTAIN SIGMA N^2 CALCULATED
        sigmaN, _ = calculate.Calculate(field_main, cg_wnoise_sum, cg_sec)
        collect_sigmaN = np.append(collect_sigmaN, sigmaN)
        
        # OBTAIN P2D AND P3D           
        ks, powspec2D, _ = get_spectrumcf(cg_wnoise_sum) 
        _, powspec2D_ave, _ = get_spectrumcf(cg_wnoise) 
        ks_n, powspec2D_noise, _ = get_spectrumcf(noiseonly) 
        powerspec3D = 2 * ks * powspec2D

        # OBTAIN RHALF AND SIGMA RHO 
        if lowpassfilter == 'no':
            brunts_sigma, _, brunts_rsq, _ = parsevals.parsevals_th(dimension, powspec2D, powerspec3D, ks, sigmaN_sq=sigmaN, targeting=target) 
        if lowpassfilter == 'yes':
            brunts_sigma, _, brunts_rsq, _= parsevals.parsevals_th(dimension, powspec2D, powerspec3D, ks, sigmaN_sq=sigmaN, the_k=ksnr) 
            
        collect_rsq = np.append(collect_rsq, brunts_rsq)
        collect_sigma = np.append(collect_sigma, brunts_sigma.item())  

        # COLLECT SPECTRAS
        # SNR
        if many == 0:
            collect_spectra = [powspec2D]
        if many > 0:
            collect_spectra = np.append(collect_spectra, [powspec2D], axis=0)  
        # SNR_AVE
        if many == 0:
            collect_spectra_ave = [powspec2D_ave]
        if many > 0:
            collect_spectra_ave = np.append(collect_spectra_ave, [powspec2D_ave], axis=0)  
        # NOISE
        if many == 0:
            collect_spectranoise = [powspec2D_noise]
        if many > 0:
            collect_spectranoise = np.append(collect_spectranoise, [powspec2D_noise], axis=0) 

    the_sigmaNmed = np.median(collect_sigmaN)
    the_rsmed = np.median(collect_rsq)
    the_sigmasmed = np.median(collect_sigma)

    the_spectrasnr = np.median(collect_spectra, axis=0)
    the_spectrasnr_ave = np.median(collect_spectra_ave, axis=0)
    the_spectranoise = np.median(collect_spectranoise, axis=0)

    return the_rsmed, the_sigmasmed, the_spectrasnr, the_spectrasnr_ave, the_spectranoise, ks, the_sigmaNmed 




# MAIN: SETUP!
# GENERAL BOX SETTINGS 
field_main = 'density'
field_two = 'cell_volume'
polytropic_index = 1 # +/-0.1, for 50mu voids

projection = True 
dimension = 1
iterations = 1000

# BOX OPTIONS: 15 total per one void size.
void_size = '50mu'
# 3 different box locations for the x los: left, center, right. and 2 more for the z los: front & back 
if 1:
    boxes = ['l_box', 'c_box', 'r_box', 'b_box', 'f_box'] 
if 0: #debug
    boxes = ['c_box'] 
# 3 different projections
if 1:  
    projection_axis_list = [0, 1, 2]  
if 0: #debug
    projection_axis_list = [1]
target_snr = [1000000, 0.01, 0.1, 1, 3, 5, 10, 100]  

# LOAD THE raw DATA
the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_512_%s_bubbles/NIF_hdf5_plt_cnt_0375'%(void_size))
# box dimensions
delta_a = 0.02  
edge_mtop = 0.3675  
edge_mbot = 0.3275
edge_mcenter = 0.3475

# IMAGES
if 0:  
    proj_2 = yt.SlicePlot(the_ds, 2, ('gas','density')) 

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

if 0:
    dx = 3.9062486030E-04
    projection_length = 2*dx  #try with 1dx 

    the_ledge = the_ds.arr([-0.15, edge_mbot, -0.15],'code_length')
    the_ledge_slice = the_ds.arr([-0.15, edge_mcenter-projection_length/2, -0.15],'code_length')
    the_redge = the_ds.arr([0.15, edge_mtop, 0.15],'code_length')
    the_redge_slice = the_ds.arr([0.15, edge_mcenter+projection_length/2, 0.15],'code_length')
    the_center =0.5*the_ledge + 0.5*the_redge
    the_region_slice = the_ds.region(the_center, the_ledge_slice, the_redge_slice)  
    the_region = the_ds.region(the_center, the_ledge, the_redge)  

    proj_1 = yt.ProjectionPlot(the_ds, 1, ('gas','density'), data_source=the_region_slice)  #slices dont' accept cgs 

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
        cfp.plot_colorbar(cmap='magma', vmin=vmin, vmax=vmax, label=label, log=True, orientation='vertical', panels=3, \
                            save='cb_location.pdf')
    stop()

if 0:
    noisetype = 'half'
    target_snr = [0.01, 0.1, 1, 3, 5, 10, 100]
    labels = ['SNR = 0.01', 'SNR = 0.1', 'SNR = 1', 'SNR = 3', 'SNR = 5', 'SNR = 10', 'SNR = 100']
    boxing = 'c_box'
    los = 1
    cg_main,_, _ = get_box(void_size, boxing, projection, los)
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
if 0:
    stop()
    print('proceed to generate/store data? c to continue')




# GETTING THE INFO WE NEED TO STORE
collect_sigmas = defaultdict(list)
collect_rsqrtsratio = defaultdict(list)
collect_rsqrt = defaultdict(list)
collect_kvals = defaultdict(list)
collect_p2dsInt = defaultdict(list)  #new
collect_sigmasNdenoised = defaultdict(list)  #new
collect_sigmasdenoised =defaultdict(list) #new
collect_rsqrttrue = defaultdict(list)
collect_spectras = defaultdict(list) 
collect_spectras_ave = defaultdict(list) 
collect_noisespectras = defaultdict(list)

# 3D data
calculate_3dsigmas = False#True
sigmas3D = []

# No noise data
no_noise = False 
if no_noise == True:
    target_snr=[1000000]
typeofnoise = 'half'
sigmaBrunt = []
sigma2D = []
rhalfs = []

# three projection axes
for eje, los in enumerate(projection_axis_list):
    # five box locations
    for num, boxing in enumerate(boxes):  
        cg_main_sum, cg_main, cg_vol = get_box(void_size, boxing, projection, los)

        if calculate_3dsigmas == True: 
            _, sigma3D = calculate.Calculate(field_main, cg_main_sum, cg_vol)
            stringing = str(sigma3D)
            numeric_part = stringing.split()[0]
            numeric_value = float(numeric_part)
            sigmas3D.append(numeric_value)
            continue

        # seven SNR
        for val,target in enumerate(target_snr):
            lowpassfilter = 'no'  
            rsqrt_med, sigma_med, collectspec, collectspec_ave, collectnoisespec, kvals2D, sigmaN_med = noise_iterations(iterations, cg_main_sum, cg_vol, lowpassfilter, target, correct=cg_main, typeofnoise=typeofnoise)
            
            sigmaBrunt.append(sigma_med)
            sigma2D.append(sigmaN_med)
            rhalfs.append(rsqrt_med)

            if no_noise == False:
                # COLLECT SPECTRAS
                if los == projection_axis_list[0] and boxing == boxes[0]:
                    collect_spectras[target] = [collectspec]
                    collect_spectras_ave[target] = [collectspec_ave]
                    collect_noisespectras[target] = [collectnoisespec]
                else: 
                    collect_spectras[target] = np.append(collect_spectras[target], [collectspec], axis=0) 
                    collect_spectras_ave[target] = np.append(collect_spectras_ave[target], [collectspec_ave], axis=0) 
                    collect_noisespectras[target] = np.append(collect_noisespectras[target], [collectnoisespec], axis=0) 

            if target != 1000000:
                rsqrt_k = []
                lowpassfilter = 'yes'  

                medspec2D = collectspec  
                medspec2D_ave = collectspec_ave  
                medspec3D = 2 * kvals2D * medspec2D
                medspec3D_ave = 2 * kvals2D * medspec2D_ave
                signalonly = collect_spectras[1000000][0]
                meannoisespec = collectnoisespec

                print('target ',target)
                #index = np.argmin(np.abs(signalonly[1:] - meannoisespec[1:]))   
                index = np.argmin(np.abs(signalonly - meannoisespec))   
                print('index ', index)
                #kcross = kvals2D[index]   
                kcross = index   
                #crossing_signalnoise = signalonly[1:][int(kcross)]
                crossing_signalnoise = signalonly[:][index]
                collect_p2dsInt[target].append(crossing_signalnoise) 
                #stop()

                if 1:
                    _, sigmaNdenoised, rsqrt_true, sigmadenoised = parsevals.parsevals_th(dimension, medspec2D, medspec3D, kvals2D, \
                                                                   sigmaN_sq=sigmaN_med, the_k=kcross, p2d_ave=medspec2D_ave, p3d_ave=medspec3D_ave)  
                rsqrt_ratio = rsqrt_true/rsqrt_med

                collect_sigmasNdenoised[target].append(sigmaNdenoised) 
                collect_sigmasdenoised[target].append(sigmadenoised)
                collect_kvals[target].append(kcross)   
                collect_rsqrttrue[target].append(rsqrt_true) 
                collect_rsqrtsratio[target].append(rsqrt_ratio) 
                collect_sigmas[target].append(sigma_med)  
                collect_rsqrt[target].append(rsqrt_med)   
                print('box snr data ready ', target,boxing,los)

# convert the default dict to a vector array 
if no_noise == False:
    data_sigmas = [*collect_sigmas.values()]
    data_rsqrts = [*collect_rsqrt.values()]
    data_kvals = [*collect_kvals.values()]
    data_rsqrttrue = [*collect_rsqrttrue.values()]
    data_rsqrtratio = [*collect_rsqrtsratio.values()]

    data_spectranoise = [*collect_spectras.values()]
    data_spectranoise_ave = [*collect_spectras_ave.values()]
    data_noisespectra = [*collect_noisespectras.values()]

    data_powintersects = [*collect_p2dsInt.values()] 
    data_sigmascolumndenoise = [*collect_sigmasNdenoised.values()] 
    data_sigmasdenoised = [*collect_sigmasdenoised.values()] 
    #stop()
    print('before storing information')
# and store the information 
if 1:
    hfivefile = '/home/552/lv9897/sflaser/h5files/15cubes15spectra_i1k_1p0_%s_512.h5'%void_size
    Fptr = h5py.File(hfivefile,'w') 
    Fptr['sigmas_%s'%void_size] = data_sigmas
    Fptr['rsqrtfull_%s'%void_size] = data_rsqrts
    Fptr['kvals_%s'%void_size] = data_kvals
    Fptr['rsqrttrue_%s'%void_size] = data_rsqrttrue
    Fptr['ratios_%s'%void_size] = data_rsqrtratio
    Fptr['kvals2D_%s'%void_size] = kvals2D

    Fptr['powintersects_%s'%void_size] = data_powintersects
    Fptr['columnddenoise_%s'%void_size] = data_sigmascolumndenoise 
    Fptr['sigmasdenoised_%s'%void_size] = data_sigmasdenoised 

    Fptr['thepowspec_%s'%void_size] = data_spectranoise 
    Fptr['thepowspec_ave_%s'%void_size] = data_spectranoise_ave 
    Fptr['thenoisepowspec_%s'%void_size] = data_noisespectra 
    Fptr.close()
    print('finished. check h5 file')
if 0:
    Fptr['sigma_3Ds_%s'%void_size] = sigmas3D
if no_noise==True:
    hfivefile = '/home/552/lv9897/sflaser/h5files/15nonoise_%s.h5'%void_size
    Fptr = h5py.File(hfivefile,'w') 
    Fptr['sigma_brunts_%s'%void_size] = sigmaBrunt
    Fptr['sigma_2Ds_%s'%void_size] = sigma2D
    Fptr['sigma_rhalfs_%s'%void_size] = rhalfs 
    Fptr.close()
    print('finished. check h5 file')
    stop()
