
from startlaser import * 
from collections import defaultdict
# README: 0 (info on data sets)
the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_50mu_bubbles/NIF_hdf5_plt_cnt_0375')

'''
to obtain SIGMA_RHO:
'''
# if field_main = any velocity, then set field_two = 'density', field_three = pressure
# README: 1 (field considered)
# SETUP:
field_main = 'density'
field_two = 'cell_volume'
field_three = None 

projection = True
dimension = 1
noise = False
noiseplot_only = False

proj_axis = 0  #or simply axis 
polytropic_index = 1 # +/-0.1, for 50um voids
the_level = 0


# CONSTRUCTING THE VOLUMES, pick:
box = 'sx_box' #and x or z 

# "SMALL" boxes considered.
delta_a = 0.0175 
edge_atop = 0.365 
edge_abot = 0.33
# "BIG" boxes ...too big..
delta_b = 0.02667  #(1.6/3)/2 
edge_btop = 0.37417 
edge_bbot = 0.32083
# "MEDIUM" boxes.
delta_m = 0.02  #(1.2/3)/2 
edge_mtop = 0.3675 
edge_mbot = 0.3275
# VARIABLE volumes
delta = 0.01


# README: 2 (in decibels)
# SIGNAL TO NOISE RATIOs & COLORS: 
target_snr = [inf]
#target_snr = [0.001, 0.01, 0.1, 1, 3, 5, 10, 100]
ksnr = [8, 9, 19, 30, 36, 44, 53]
#target_colors = ['orange','y','magenta','blue','green','grey'] #'black']
target_colors = ['blue','magenta','green']

iterations = 1
noiseat = 'none'  # when to add noise, at 'image' or at 'dft - no noise really...'
lowpassfilter = 'no'  # yes:low k values, no:all k values

method = 0  # method to obtain sigma N - 0: via PDF, 1: via np.std(), 2: rms 
imaging = 'no'  # make images?
width = 2*delta_m  # for images


collect = defaultdict(list)
# FOR PLOTS 
if 1:
    fig,ax = plt.subplots(1,1)
    fig.subplots_adjust(wspace=0, hspace=0)

for pull in range(6): 
    pull = pull + 1
    if box == 'sx_box':
        the_ledge = the_ds.arr([-delta*pull, edge_mbot, -delta_m],'code_length')
        the_redge = the_ds.arr([delta*pull, edge_mtop, delta_m],'code_length')
    if box == 'sz_box':
        the_ledge = the_ds.arr([-delta_m, edge_mbot, -delta*pull],'code_length')
        the_redge = the_ds.arr([delta_m, edge_mtop, delta*pull],'code_length')


    # THE PIPELINE FOR THE CHOICES OF VOLUME
    the_center =0.5*the_ledge + 0.5*the_redge
    mindx = the_ds.index.get_smallest_dx()
    size = ((the_redge - the_ledge)/mindx).astype('int') 

    the_region = the_ds.region(the_center, the_ledge, the_redge)

    section = extract_s1.extract_section_1(the_ds,the_level,the_ledge,size)
    section_blur = motion_blur.Motion_Blur()  #tbd; simply returns
    cg_first,cg_second,cg_third = project.Project(section,proj_axis,projection,field_main,field_two,field_three)


    for val,target in enumerate(target_snr):
        for many in range(iterations):  

            # ADD NOISE TO THE IMAGE
            if 1:  
                if noiseat == 'image':
                    cg_wnoise, noiseonly = noiseimage.add_noise(cg_first, target)
                if noiseat == 'none':  
                    cg_wnoise = cg_first
                    noiseonly = 0
            
            # PLOT POWERSPECTRUMS
            if 0:  
                k1, spec1, _, _ = powerspec.powerspec_lowpass(cg_wnoise)
                #k2, spec2, _, _ = powerspec.powerspectrum_ii(cg_first)
                k3, spec3, _, _ = powerspec.powerspec_lowpass(noiseonly, justnoise=noise)

                ax.plot(k1, spec1, c=target_colors[val], alpha=0.5)
                #ax.plot(k2, spec2, c='b')
                ax.plot(k3, spec3, c='k', alpha=0.5)
                ax.set(yscale='log')#, xscale='log')

            # IMAGE THE SIGNAL, NOISE, SIGNAL PLUS NOISE
            if 0:  
                cfp.plot_map(cg_first.v,save='signal')
                cfp.plot_map(cg_wnoise,save='signal_noise.png')
                cfp.plot_map(noiseonly,save='noise.png')
                stop()

            # OBTAIN SIGMA N^2
            if 1:  #calculated
                sigma_sq = calculate.Calculate(field_main, the_region, cg_wnoise, cg_second) 
            if 0:  #from PDF
                bins, cbins, pdf = make_pdfNfit.Make_PDFnFIT(field_main, projection, proj_axis, cg_wnoise ,cg_third)  #EDIT; split into definitions 
                sigma_sq = make_sigmasquared.Make_SigmaSquared(field_main, bins, cbins, pdf, cg_wnoise, method)
            if projection == False: 
                sigma_sq = None

            # OBTAIN P2D AND P3D           
            if 1:   
                if lowpassfilter == 'no':
                    four_amps, powerspec2D, powerspec3D, infxyz, ks = powerspec.powerspec(cg_wnoise, dimension, box, target, many)    
                if lowpassfilter == 'yes':
                    ks, powerspec2D, powerspec3D, infxyz = powerspec.powerspec_lowpass(cg_wnoise, the_k=ksnr[val], lowpass_filter=lowpassfilter)
            if 0:  # if adding noise to DFT: DELETE!
                noisy_powerspec2D, noisy_powerspec3D, noise_sign, noise_signs = signaltonoise.signaltonoise(cg_first, dimension, target, four_amps, powerspec2D, powerspec3D) 
                # README 3. (thoughts on making pdfs)
            
            # OBTAIN SIGMA RHO
            if 1:
                sigma_rho, brunts_rsq = parsevals.parsevals_th(dimension, powerspec2D, powerspec3D, infxyz, ks, sigma_sq)
                #collect[target].append(brunts_rsq)   #plotting brunts_rsq  
                collect[target].append(sigma_rho)    #plotting sigma_rho 

    # SAVE 
    if 0:  #...the powerspectrum plots
        plt.ylim(1e0,1e7)
        plt.xlim(40,75)
        fig.savefig('pd2_snr_noise_logy_cbox_510100_zoom_.png')#%target)


# PLOTS & FITS: for SX and SZ boxing.
if 1:
    if box == 'sx_box':
        the_truey = [0.31001771892983476, 0.3078297687398883, 0.30420254190841384, 0.2995425654416588, 0.29491911865113646, 0.2943182444076409]
        the_trueylog = np.log10(the_truey) 
    if box == 'sz_box':
        the_truey = [0.29975269098689805, 0.3078297687398883, 0.318353770910557, 0.3112775959416291, 0.30599288334855623, 0.3004886389624372]
        the_trueylog = np.log10(the_truey) 

    the_d = np.linspace(0.02, 0.12, pull, endpoint=True)  
    the_l = delta_m * 2

    the_x = the_d/the_l # aspect ratio 
    the_xlog = np.log10(the_x)

    the_prey = [l for l in collect.values()]
    emptylist = []
    the_y = np.append(emptylist,the_prey[0])

    the_ylog = np.log10(the_y) 
    the_yerr = np.std(the_ylog)
 
    y_errlog= 0.1 * the_ylog
    y_err = 0.1 * the_y
    y_err_arr = [y_err,y_err]


    # POSSIBLE FITTINGS
    def afunct(x, a, b):
        #y = a * x ** b
        y = a * x + b
        return y
    def bfunct(x, a, b):
        y = 1.5 * x ** b  #this may have thought to be a 1.5...see slack, and it was really -0.5
        return y
    if 0:
        param, param_cov = curve_fit(bfunct, the_x, the_y)
        plt.plot(the_x, bfunct(the_x, *param), 'b--')

    # FITTING A CURVE
    if 1:
        rets = cfp.fit(afunct, the_xlog, the_ylog, yerr=the_yerr)
        
        the_xrecipe = np.linspace(the_xlog.min(),the_xlog.max(),num=len(the_xlog))
        the_xten = 10**the_xrecipe
        the_yrecipe = (rets.popt[0]*the_xrecipe + rets.popt[1]) #afunct(the_xlog, *rets.popt)
        the_yten = 10**the_yrecipe

        plt.plot(the_xten, the_yten, 'g--')
        #cfp.plot(the_xten, the_yten, color='g', linestyle='dashed')
        print('param cf vals ',rets.popt)

    # CONTINUE WITH CFPACK 
    if 1:
        cfp.plot(x=the_x, y=the_truey, type='scatter', color='orange')
        cfp.plot(x=the_x, y=the_y, yerr=y_err_arr, marker = 'o', xlabel='aspect ratio', ylabel=r'$\sigma_{\rho/\rho_o}$', xlog=True, ylog=True, \
                save='sigma_rho_%s_aspect_logaxes_fit_errbarsnotlogged_cfp.png'%(box))
    # OR TRADITIONALLY...
    if 0:
        ax.errorbar(the_x, the_y, yerr=0.1)  #plt.scatter(the_x, the_y, color='b')
        plt.scatter(the_x, the_truey, color='orange')
        ax.set_yscale('log')
        ax.set_xscale('log')
        plt.xlabel('aspect ratio')
        plt.ylabel(r'$\sigma_{\rho/\rho_o}$')
        plt.savefig('sigma_rho_%s_aspect_logaxes_fit_errbarstest.png'%(box))
        plt.clf() 




# CALLING BOXPLOTS; README 4 (an example). README 5 (plotting, dealing with BIG numbers in axis).
# TRUE values - l_c_r: 0.269, 0.308, 0.265; b_c_f: 0.288, 0.308(yay), 0.243 
if 0:  
    labels, data = kcollect.keys(), kcollect.values()  #PICK THE RIGHT DICTIONARY
    plt.boxplot(data, showfliers=True, showmeans=True)

    the_mean = 0.303531
    plt.plot([1,len(kmatch)],[the_mean,the_mean],'--',c='k',alpha=0.7)
    plt.plot([1,len(kmatch)],[the_rsallmean,the_rsallmean],'--',c=[0.5]*4)
    '''
    if box == 'l_box':
        plt.plot([1,6],[0.269,0.269],'--',c=[0.5]*4)
    if box == 'c_box':
        plt.plot([1,6],[0.308,0.308],'--',c=[0.5]*4)
    if box == 'r_box':
        plt.plot([1,6],[0.265,0.265],'--',c=[0.5]*4)
    '''
    plt.xlabel('k')
    #plt.ylabel(r'$\sigma_{\rho/\rho_o} $')
    #plt.ylabel(r'$R^{1/2}$')
    plt.xticks(range(1, len(labels) + 1), labels)
    #plt.yscale('log')
    #plt.ylim(0.10,0.39)
    #plt.savefig('boxplot_1000n_snr_%s_xpos_image_bruntsrsqVSsnr_001_fixed'%box)
    plt.savefig('best_kvalue_%s_00%d_winf_wresp'%(box,target))
    print(box)
    plt.clf()

