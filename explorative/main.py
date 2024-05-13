
from startlaser import * 
from collections import defaultdict
# README: 0 (info on data sets)

void_size = '50um'
if void_size=='12p5um': 
    the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_12p5mu_bubbles/NIF_hdf5_plt_cnt_0375')
if void_size=='50um': 
    the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_50mu_bubbles/NIF_hdf5_plt_cnt_0375')
if void_size=='100um':
    the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_100mu_bubbles/NIF_hdf5_plt_cnt_0375')

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
noise = True
noiseplot_only = False

proj_axis = 0  #or simply axis 
polytropic_index = 1 # +/-0.1, for 50um voids
the_level = 0

# CONSTRUCTING THE VOLUMES, pick:
boxes = ['l_box','c_box','r_box'] #['b_box','c_box','f_box'] 

# "SMALL" boxes considered.
delta_a = 0.0175 
edge_atop = 0.365 
edge_abot = 0.33
# "BIG" boxes ...too big..
delta_b = 0.02667  #(1.6/3)/2 
edge_btop = 0.37417 
edge_bbot = 0.32083

# "MEDIUM" boxes & aspect test.
delta_a = 0.02  #(1.2/3)/2
delta_m = 0.02  #0.03 for aspects  
edge_mtop = 0.3675 
edge_mbot = 0.3275


# README: 2 (in decibels)
# SIGNAL TO NOISE RATIOs & COLORS: :w

#target_snr = [inf]
target_snr = [0.01, 0.1, 1, 3, 5, 10, 100]
ksnr = [6, 8, 18, 29, 35, 44, 52]  #right: current, lef:[5, 7, 16, 25, 30, 38, 52], center: ksnr = [8, 8, 9, 19, 30, 36, 44, 52]
target_colors = ['orange','y','magenta','blue','green','cyan','red'] #'black']
#target_colors = ['blue','magenta','green']

iterations = 100
noiseat = 'image'  # when to add noise, at 'image' or at 'dft - no noise really...'
lowpassfilter = 'no'  # yes:low k values, no:all k values

method = 0  # method to obtain sigma N - 0: via PDF, 1: via np.std(), 2: rms 
imaging = 'no'  # make images?
width = 2*delta_m  # for images


# OVERLAYING 3 IND BOXES RESULTS
if 0:
    fig,ax = plt.subplots(1,1)
    fig.subplots_adjust(wspace=0, hspace=0)
    collect = defaultdict(list)

# THREE BOXES
limit = np.zeros(4) #for README 5
for num, box in enumerate(boxes): 
    if box == 'c_box':
        the_ledge = the_ds.arr([-delta_a, edge_mbot, -delta_a],'code_length')
        the_redge = the_ds.arr([delta_a, edge_mtop, delta_a],'code_length')

    if box == 'r_box':
        the_ledge = the_ds.arr([delta_m, edge_mbot, -delta_m],'code_length')
        the_redge = the_ds.arr([delta_m+2*delta_m, edge_mtop, delta_m],'code_length')
        #the_ledge = the_ds.arr([0, edge_mbot, -delta_a],'code_length')
        #the_redge = the_ds.arr([delta_m+2*delta_m, edge_mtop, delta_a],'code_length')
    if box == 'f_box':
        the_ledge = the_ds.arr([-delta_m, edge_mbot, delta_m],'code_length')
        the_redge = the_ds.arr([delta_m, edge_mtop, delta_m+2*delta_m],'code_length')
        #the_ledge = the_ds.arr([-delta_a, edge_mbot, 0],'code_length')
        #the_redge = the_ds.arr([delta_a, edge_mtop, delta_m+2*delta_m],'code_length')

    if box == 'l_box':
        the_ledge = the_ds.arr([-delta_m-2*delta_m, edge_mbot, -delta_m],'code_length')
        the_redge = the_ds.arr([-delta_m, edge_mtop, delta_m],'code_length')
        #the_ledge = the_ds.arr([-delta_m-2*delta_m, edge_mbot, -delta_a],'code_length')
        #the_redge = the_ds.arr([0, edge_mtop, delta_a],'code_length')
    if box == 'b_box':
        the_ledge = the_ds.arr([-delta_m, edge_mbot, -delta_m-2*delta_m],'code_length')
        the_redge = the_ds.arr([delta_m, edge_mtop, -delta_m],'code_length')
        #the_ledge = the_ds.arr([-delta_a, edge_mbot, -delta_m-2*delta_m],'code_length')
        #the_redge = the_ds.arr([delta_a, edge_mtop, 0],'code_length')
    '''
    if box == 'notabox':
        # original analysis region:  #1.2 mm × 0.34 mm × 1.2 mm.
        the_ledge = the_ds.arr([-0.06,0.33,-0.06],'code_length')
        the_redge = the_ds.arr([0.06,0.365,0.06],'code_length')
    '''

    if 1:  
        if 1:
            collect = defaultdict(list)
        # FOR PLOTS 
        if 1:
            fig,ax = plt.subplots(1,1)
            fig.subplots_adjust(wspace=0, hspace=0)

        # THE PIPELINE FOR ALL CHOICES OF VOLUME
        the_center =0.5*the_ledge + 0.5*the_redge
        mindx = the_ds.index.get_smallest_dx()
        size = ((the_redge - the_ledge)/mindx).astype('int') 

        the_region = the_ds.region(the_center, the_ledge, the_redge)

        section = extract_s1.extract_section_1(the_ds,the_level,the_ledge,size)
        section_blur = motion_blur.Motion_Blur()  #tbd; simply returns
        cg_first,cg_second,cg_third = project.Project(section,proj_axis,projection,field_main,field_two,field_three)

        if 0:  # EARLY IMAGES: EDIT 
            images = make_images.Make_Images(the_ds,cg_first,proj_axis,box,the_center,width) 
        if 0:
            mom_map = momentmap.moment_map(cg_first)  #a brief exploration, maybe discard

        
        rsqrt = []
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
                if 1:  
                    k1, spec1, _, _ = powerspec.powerspec_lowpass(cg_wnoise)
                    #k2, spec2, _, _ = powerspec.powerspectrum_ii(cg_first)
                    k3, spec3, _, _ = powerspec.powerspec_lowpass(noiseonly, justnoise=noise)

                    ax.plot(k1, spec1, c=target_colors[val], alpha=0.5)
                    #ax.plot(k2, spec2, c='b')
                    ax.plot(k3, spec3, c='k', alpha=0.5)
                    ax.set(yscale='log')#, xscale='log')

                # IMAGE THE SIGNAL, NOISE, SIGNAL PLUS NOISE
                if 0:  
                    cfp.plot_map(cg_first.v,save='signal_%s'%box)#, vmin=20, vmax=49)
                    #cfp.plot_map(cg_wnoise,save='signal_noise.png')
                    #cfp.plot_map(noiseonly,save='noise.png')
                    #stop()

                # OBTAIN SIGMA N^2
                if 0:  #calculated
                    sigma_sq = calculate.Calculate(field_main, the_region, cg_wnoise, cg_second) 
                if 0:  #from PDF
                    bins, cbins, pdf = make_pdfNfit.Make_PDFnFIT(field_main, projection, proj_axis, cg_wnoise ,cg_third)  #EDIT; split into definitions 
                    sigma_sq = make_sigmasquared.Make_SigmaSquared(field_main, bins, cbins, pdf, cg_wnoise, method)
                if projection == False: 
                    sigma_sq = None

                # OBTAIN P2D AND P3D           
                if 0:   
                    if lowpassfilter == 'no':
                        four_amps, powerspec2D, powerspec3D, infxyz, ks = powerspec.powerspec(cg_wnoise, dimension, box, target, many)    
                    if lowpassfilter == 'yes':
                        ks, powerspec2D, powerspec3D, infxyz = powerspec.powerspec_lowpass(cg_wnoise, the_k=ksnr[val], lowpass_filter=lowpassfilter)
                if 0:  # if adding noise to DFT: DELETE!
                    noisy_powerspec2D, noisy_powerspec3D, noise_sign, noise_signs = signaltonoise.signaltonoise(cg_first, dimension, target, four_amps, powerspec2D, powerspec3D) 
                    # README 3. (thoughts on making pdfs)
                
                # OBTAIN SIGMA RHO
                if 0:
                    sigma_rho, sigma_rho2d, brunts_rsq = parsevals.parsevals_th(dimension, powerspec2D, powerspec3D, infxyz, ks, sigma_sq)
                    if target == inf:
                        print('at the original signal')
                    #collect[box].append(sigma_rho)   #plotting brunts_rsq  
                    collect[target].append(brunts_rsq)   #plotting brunts_rsq  
                    #collect[target].append(sigma_rho)    #plotting sigma_rho 


                
            # LOW PASS FILTERING & OTHERS.. 
            if 0:
                the_rs = [l for l in collect.values()]
                emptylist = []
                the_rsall = np.append(emptylist,the_rs[val])

                # for the cfpack fit
                the_rsallmed = np.median(the_rsall)

                the_rsalllogmed = np.log10(the_rsallmed)
                the_rsalllog = np.log10(the_rsall)
                the_rsallmedlog = np.median(the_rsalllog) 
                the_rsalllogerr = np.std(the_rsalllog)

                rsqrt = np.append(rsqrt, the_rsallmed)
                #stop() #logmed ~ medlog, check [ ]

    # SAVE THE POWERSPECTRUM PLOTS 
    if 1:  
        plt.ylim(1e1,1e11)
        #plt.xlim(40,75)
        fig.savefig('pd2_snr_sn_wn_logy_%s_%d.png'%(box,proj_axis))
    print('okay')

    # to learn which are the crossing k values - then these get plotted at ''CALLING BOXPLOTS'' 
    if 0:
        kcollect = defaultdict(list)  #CAREFUL?
        kmatch = np.linspace(1,16,16)  #target 0.1
        #kmatch = np.linspace(10,25,16)  #target 1
        #kmatch = np.linspace(20,35,16)  #target:5,3
        #kmatch = np.linspace(34,54,11)  #target 10
        #kmatch = np.linspace(50,60,11)   #target 100
        target = 0.01
        iterations = 1000
        for val,ksnr in enumerate(kmatch):
            for many in range(iterations):  
                cg_wnoise, noiseonly = noiseimage.add_noise(cg_first, target)
                sigma_sq = calculate.Calculate(field_main, the_region, cg_wnoise, cg_second) 
                ks, powerspec2D, powerspec3D, infxyz = powerspec.powerspec_lowpass(cg_wnoise, the_k=ksnr, lowpass_filter='yes')
                sigma_rho, brunts_rsq = parsevals.parsevals_th(dimension, powerspec2D, powerspec3D, infxyz, ks, sigma_sq)
                kcollect[ksnr].append(brunts_rsq)   

    # crossing k values vs SNR -raw
    if 0:
        ydata = [8, 8, 9, 19, 29, 36, 44, 53]
        ydatalog = np.log10(ydata)
        xdata = target_snr
        #xdatalog = np.log10(xdata)
        if 0:
            param, param_cov = curve_fit(cfunct, xdata, ydata)
            plt.plot(xdata, cfunct(xdata, *param), 'b--')
        plt.scatter(xdata, ydata)
        ax.set(xscale='log', yscale='log')
        plt.ylabel('Crossing k values')
        plt.xlabel('SNR')
        plt.savefig('crossing_kvalues_vs_SNR_%s'%(box))
        #stop()
        plt.clf() 

    # POSSIBLE FIT FUNCTIONS
    def cfunct(x, a, b, c, d):
        y = a * np.tanh((c*x)+d) + b
        return y

    def dfunct(x, a, b):
        y = 1/(1 + a*np.exp(-x+b)) 
        return y

    # r^{1/2}true / r^{1/2}full vs SNR - raw
    rhalf_ctrue = 0.3035  #0.3035309747962126
    rhalf_ltrue = 0.3697  #0.3697364679188141
    rhalf_rtrue = 0.3385  #0.33851156031827845
    rhalf_btrue = 0.3274  #0.3273638724473809
    rhalf_ftrue = 0.3281  #0.3280843502158977
    if 0:

        ytrue_cbox = [0.30200493, 0.3006561,  0.30287845, 0.30259622, 0.30320516, 0.30334879, 0.30352811]
        ytrue_lbox = [0.36824365, 0.3592389,  0.36746576, 0.36903979, 0.37024448, 0.37046739, 0.3696961] 
        ytrue_rbox = [0.34350638, 0.3291993,  0.33969969, 0.33873154, 0.33918174, 0.33809821, 0.3384887]

        ydata = rsqrt/ytrue_cbox
        ydatalog = np.log10(ydata)
        xdata = target_snr
        xdatalog = np.log10(xdata)
        if 1:  #tanh
            parametros = defaultdict(list)
            #parametros['d'].append(0)
            #parametros['d'].append(2)
            #parametros['d'].append(4)
            rets = cfp.fit(cfunct, xdatalog, ydatalog)
            #stop()
            the_xrecipe = np.linspace(xdatalog.min(),xdatalog.max(),num=500)  # more numbers for the x fitting! ,was len(xdatalog)
            the_xten = 10**the_xrecipe
            the_yrecipe = cfunct(the_xrecipe, *rets.popt)
            the_yten = 10**the_yrecipe

            plt.plot(the_xten, the_yten, 'g--')
            print('param cf vals ',rets.popt) 
        if 0:  #sigmoid
            rets = cfp.fit(dfunct, xdatalog, ydatalog)
            the_xrecipe = np.linspace(xdatalog.min(),xdatalog.max(),num=len(xdatalog))
            the_xten = 10**the_xrecipe
            the_yrecipe = (1 / (1+ rets.popt[0] * np.exp(-1*the_xrecipe)))
            the_yten = 10**the_yrecipe
            plt.plot(the_xten, the_yten, 'g--')
            print('param cf vals ',rets.popt) 

        rhalf_true = rhalf_ltrue 
        plt.scatter(xdata, ydata)  
        #plt.plot([1e-2,100],[rhalf_true,rhalf_true],'--',c='grey')
        ax.set(xscale='log',yscale='log')
        plt.ylabel(r'$R^{1/2}_{true}/R^{1/2}_{full}$') #plt.ylabel(r'$R^{1/2}_{full}$')
        plt.xlabel('SNR')
        #ax.set(ylim=[2e-1,1.1e0])  #ax.set(xlim=[1e-3,2e2], ylim=[1e-1,4e-1])
        plt.savefig('rsqrt_ratio_vs_SNR_tanhx_multplus_%s'%(box))  #rhalf_full, rhalf_true, rhalf_ratio
        plt.clf()



    # CALLING BOXPLOTS; README 4 (an example). README 5 (plotting, dealing with BIG numbers in axis).
    # TRUE values - l_c_r: 0.269, 0.308, 0.265; b_c_f: 0.288, 0.308(yay), 0.243 
    if 0:  
        labels, data = collect.keys(), collect.values()  #PICK THE RIGHT DICTIONARY
        plt.boxplot(data, showfliers=True, showmeans=False)

        the_mean = rhalf_ctrue  #rsqrt at inf = pure signal 
        plt.plot([1,len(labels)],[the_mean,the_mean],'--',c='k',alpha=0.7)
        #plt.plot([1,len(kmatch)],[the_rsallmean,the_rsallmean],'--',c=[0.5]*4)
        '''
        if box == 'l_box':
            plt.plot([1,6],[0.269,0.269],'--',c=[0.5]*4)
        if box == 'c_box':
            plt.plot([1,6],[0.308,0.308],'--',c=[0.5]*4)
        if box == 'r_box':
            plt.plot([1,6],[0.265,0.265],'--',c=[0.5]*4)
        '''
        plt.xlabel('SNR')
        #plt.ylabel(r'$\sigma_{\rho/\rho_o} $')
        plt.ylabel(r'$R^{1/2}$')
        plt.xticks(range(1, len(labels) + 1), labels)
        #plt.yscale('log')
        #plt.ylim(0.10,0.39)
        #plt.savefig('boxplot_1000n_snr_%s_xpos_image_bruntsrsqVSsnr_001_fixed'%box)
        #plt.savefig('best_kvalue_%s_00%d_winf'%(box,target))
        plt.savefig('rsqrt_iterations_%s_winf'%(box))
        print(box)
        plt.clf()





# INDEPENDENT BOXES FOR ASPECT RATIOS
if 0:
    the_rs = [l for l in collect.values()]
    emptylist = []
    the_sigmas = np.append(emptylist,the_rs)
    the_3DzSigmas_cube = [0.2880937570662445,0.30782976873988827, 0.24301533336823522]
    the_3DxSigmas_cube = [0.269150311726313, 0.30782976873988827, 0.2654539963372979]
    the_3DzSigmas_aspect = [0.30397296809676416, 0.318353770910557, 0.2807981869925865]  #[0.246715430918834, 0.318353770910557, 0.23660419455366774]
    the_3DxSigmas_aspect = [0.2895337967904341, 0.30420254190841384, 0.326483876406291]
    the_d_cube = np.linspace(-0.04, 0.04, 3, endpoint=True)  
    the_d = np.linspace(-0.03, 0.03, 3, endpoint=True)  
    y_err_up = 0.1*the_sigmas
    y_err_low = 0.1*the_sigmas 
    y_err_arr = [y_err_low,y_err_up]
    if 1:
        ax.scatter(the_d, the_sigmas, alpha=0.5)
        #ax.errorbar(the_d, the_sigmas, yerr=y_err_arr, fmt='o')
        ax.scatter(the_d, the_3DzSigmas_aspect, c='orange', alpha=0.5)
        plt.ylim(0.2207,0.4855) #Z
        #plt.ylim(0.2250,0.4780) #X 0.2295,0.4780
        #plt.yscale('log')
        plt.ylabel(r'$\sigma_{\rho/\rho_o}$')
        plt.xlabel('the volume centers')
        fig.savefig('Z_aspect_retry.png')

'''
to obtain SIGMA_V:
'''
'''
field = 'velocity'
proj_axis_vel = 2
#section_v = extract_s2.extract_section_2(the_ds2, the_level,the_ledge2,size2)  #for the more suitable region; using section for now
project_vel = project.Project(section,field,proj_axis_vel,the_center,the_redge,the_ledge)
#ProbablyOtherPhysicsTBD
bins_v, cbins_v, pdf_v = make_pdf.Make_PDF(project_vel,field) 
sigma_v = make_sigmavel.Make_SigmaVel(field_main,bins_v,cbins_v,pdf_v) 
'''

'''
TO OBTAIN b: 
'''
if 0: 
    obtain_b.Obtain_b()  #sigma_rho,sigma_v,polytropic_index

