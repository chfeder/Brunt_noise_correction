
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
noise = True
noiseplot_only = False

proj_axis = 1  #or simply axis 
polytropic_index = 1 # +/-0.1, for 50um voids
the_level = 0

# CONSTRUCTING THE VOLUMES, pick:
boxes = ['c_box'] #['b_box','c_box','f_box'] #['l_box','c_box','r_box'] ; 

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
#target_snr = [inf]
target_snr = [0.001, 0.01, 0.1, 1, 3, 5, 10, 100]
ksnr = [8, 8, 9, 19, 30, 36, 44, 53]
#target_colors = ['orange','y','magenta','blue','green','grey'] #'black']
target_colors = ['blue','magenta','green']

iterations = 1000
noiseat = 'image'  # when to add noise, at 'image' or at 'dft - no noise really...'
lowpassfilter = 'yes'  # yes:low k values, no:all k values

method = 0  # method to obtain sigma N - 0: via PDF, 1: via np.std(), 2: rms 
imaging = 'no'  # make images?
width = 2*delta_m  # for images


# THREE BOXES
limit = np.zeros(4) #for README 5
for num, box in enumerate(boxes): 
    if box == 'c_box':
        the_ledge = the_ds.arr([-delta_m, edge_mbot, -delta_m],'code_length')
        the_redge = the_ds.arr([delta_m, edge_mtop, delta_m],'code_length')

    if box == 'r_box':
        the_ledge = the_ds.arr([delta_m, edge_mbot, -delta_m],'code_length')
        the_redge = the_ds.arr([delta_m+2*delta_m, edge_mtop, delta_m],'code_length')
    if box == 'f_box':
        the_ledge = the_ds.arr([-delta_m, edge_mbot, delta_m],'code_length')
        the_redge = the_ds.arr([delta_m, edge_mtop, delta_m+2*delta_m],'code_length')

    if box == 'l_box':
        the_ledge = the_ds.arr([-delta_m-2*delta_m, edge_mbot, -delta_m],'code_length')
        the_redge = the_ds.arr([-delta_m, edge_mtop, delta_m],'code_length')
    if box == 'b_box':
        the_ledge = the_ds.arr([-delta_m, edge_mbot, -delta_m-2*delta_m],'code_length')
        the_redge = the_ds.arr([delta_m, edge_mtop, -delta_m],'code_length')
    '''
    if box == 'notabox':
        # original analysis region:  #1.2 mm × 0.34 mm × 1.2 mm.
        the_ledge = the_ds.arr([-0.06,0.33,-0.06],'code_length')
        the_redge = the_ds.arr([0.06,0.365,0.06],'code_length')
    '''

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


        rhalf = []
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
                    collect[target].append(brunts_rsq)   #plotting brunts_rsq  
                    #collect[target].append(sigma_rho)    #plotting sigma_rho 


                
            # LOW PASS FILTERING & OTHERS.. 
            if 1:
                the_rs = [l for l in collect.values()]
                emptylist = []
                the_rsall = np.append(emptylist,the_rs[val])

                # for the cfpack fit
                the_rsallmed = np.median(the_rsall)

                the_rsalllogmed = np.log10(the_rsallmed)
                the_rsalllog = np.log10(the_rsall)
                the_rsallmedlog = np.median(the_rsalllog) 
                the_rsalllogerr = np.std(the_rsalllog)

                rhalf = np.append(rhalf, the_rsallmed)
                #stop() #logmed ~ medlog, check [ ]

    # SAVE THE POWERSPECTRUM PLOTS 
    if 0:  
        plt.ylim(1e0,1e7)
        plt.xlim(40,75)
        fig.savefig('pd2_snr_noise_logy_cbox_510100_zoom_.png')#%target)
    print('okay')
    
    if 0:
        kcollect = defaultdict(list)  #CAREFUL?
        kmatch = np.linspace(5,20,16)  #target 0.1
        #kmatch = np.linspace(15,30,16)  #target 1
        #kmatch = np.linspace(25,40,16)  #target 3
        #kmatch = np.linspace(34,54,11)  #target 5,10
        #kmatch = np.linspace(50,60,11)
        target = 0.001 
        iterations = 1000
        for val,ksnr in enumerate(kmatch):
            for many in range(iterations):  
                cg_wnoise, noiseonly = noiseimage.add_noise(cg_first, target)
                sigma_sq = calculate.Calculate(field_main, the_region, cg_wnoise, cg_second) 
                ks, powerspec2D, powerspec3D, infxyz = powerspec.powerspec_lowpass(cg_wnoise, the_k=ksnr, lowpass_filter='yes')
                sigma_rho, brunts_rsq = parsevals.parsevals_th(dimension, powerspec2D, powerspec3D, infxyz, ks, sigma_sq)
                kcollect[ksnr].append(brunts_rsq)   


    def cfunct(x, a, b):
        y = a * x**(1/b)
        return y
    # crossing k values vs SNR -raw
    if 0:
        ydata = [8, 8, 9, 19, 29, 36, 44, 53]
        ydatalog = np.log10(ydata)
        xdata = target_snr
        #xdatalog = np.log10(xdata)
        if 0:
            param, param_cov = curve_fit(cfunct, xdata, ydata)
            plt.plot(xdata, cfunct(xdata, *param), 'b--')
        plt.scatter(xdata[1:], ydata[1:])
        ax.set(xscale='log', yscale='log')
        plt.ylabel('Crossing k values')
        plt.xlabel('SNR')
        plt.savefig('crossing_kvalues_vs_SNR_%s'%(box))
        #stop()
        plt.clf() 
    # r^{1/2}true / r^{1/2}full vs SNR - raw
    if 1:
        ydata = rhalf  
        #ydata = [0.3731, 0.3734, 0.3747, 0.4950, 0.7874, 0.9009, 0.9718, 1]
        ydatalog = np.log10(ydata)
        xdata = target_snr
        xdatalog = np.log10(xdata)
        #stop()
        if 0:  # WILL BE EDITED
            param, param_cov = curve_fit(cfunct, xdata, ydata)
            plt.plot(xdata, cfunct(xdata, *param), 'b--')

        plt.scatter(xdata[1:], ydata[1:])  #first value throws an error
        ax.set(xscale='log', yscale='log')
        plt.ylabel(r'$R^{1/2}_{true}$') #plt.ylabel(r'$R^{1/2}_{true} / R^{1/2}_{full}$') #plt.ylabel(r'$R^{1/2}_{full}$')
        plt.xlabel('SNR')
        ax.set(xscale='log', yscale='log')
        #ax.set(xlim=[1e-3,2e2], ylim=[1e-1,4e-1])
        plt.savefig('rhalf_true_vs_SNR_%s'%(box))  #rhalf_full, rhalf_true, rhalf_ratio
        #stop()
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

