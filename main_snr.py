
from startlaser import * 
from collections import defaultdict


# ALL THE OPTIONS
# 3 different void sizes
void_size = '50um'
if void_size=='12p5um': 
    the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_12p5mu_bubbles/NIF_hdf5_plt_cnt_0375')
if void_size=='50um': 
    the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_50mu_bubbles/NIF_hdf5_plt_cnt_0375')
if void_size=='100um':
    the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_100mu_bubbles/NIF_hdf5_plt_cnt_0375')

# 3 different boxes for the x los: could combine these depending on final presentation 
if 1:
    boxes_list = [['l_box'],['c_box'],['r_box']] 
    boxes = boxes_list[1]
# 3 different boxes for the z los 
if 0:
    boxes = ['b_box','c_box','f_box'] 
# notes: could fit 4 more in the diagonals... do slices to check [  ]

# 3 different projections
proj_axis_list = [0,1,2] 
proj_axis = proj_axis_list[0] 


field_main = 'density'
field_two = 'cell_volume'
field_three = None 
polytropic_index = 1 # +/-0.1, for 50um voids

projection = True 
dimension = 1
noise = True
iterations = 1000
noiseat = 'image'  
lowpassfilter = 'no'  # yes:low k values, no:all k values

# SIGNAL TO NOISE RATIOs & COLORS: :w
if 0:  #pure signal
    target_snr = [inf]
if 1:  #noise levels
    target_snr = [0.01, 0.1, 1, 3, 5, 10, 100, inf]
target_colors = ['orange','y','magenta','blue','green','cyan','red'] 

# k values to recover rsqrt_true, these were picked by hand, automate? [ ]    
ksnr = [6, 8, 18, 29, 35, 44, 52]  #right, left:[5, 7, 16, 25, 30, 38, 52], center: ksnr = [8, 8, 9, 19, 30, 36, 44, 52]


the_level = 0
# for the 'medium' boxes 
delta_a = 0.02  #(1.2/3)/2
delta_m = delta_a  
edge_mtop = 0.3675 
edge_mbot = 0.3275

# re-consider the loop! over which options
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

    # THE PIPELINE FOR ALL CHOICES OF VOLUME
    the_center =0.5*the_ledge + 0.5*the_redge
    mindx = the_ds.index.get_smallest_dx()
    size = ((the_redge - the_ledge)/mindx).astype('int') 

    the_region = the_ds.region(the_center, the_ledge, the_redge)

    section = extract_s1.extract_section_1(the_ds,the_level,the_ledge,size)
    section_blur = motion_blur.Motion_Blur()  #tbd; simply returns
    cg_first,cg_second,cg_third = project.Project(section,proj_axis,projection,field_main,field_two,field_three)

    if 0:  
        if 1:
            collect = defaultdict(list)
        # FOR PLOTS 
        if 1:
            fig,ax = plt.subplots(1,1)
            fig.subplots_adjust(wspace=0, hspace=0)


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
                    cfp.plot_map(cg_first.v,save='signal_%s'%box)#, vmin=20, vmax=49)
                    #cfp.plot_map(cg_wnoise,save='signal_noise.png')
                    #cfp.plot_map(noiseonly,save='noise.png')
                    #stop()

                # OBTAIN SIGMA N^2
                if 1:  #calculated
                    sigma_sq = calculate.Calculate(field_main, the_region, cg_wnoise, cg_second) 
                if 0:  #from PDF
                    bins, cbins, pdf = make_pdfNfit.Make_PDFnFIT(field_main, projection, proj_axis, cg_wnoise ,cg_third)  #EDIT; split into definitions 
                    method = 0
                    sigma_sq = make_sigmasquared.Make_SigmaSquared(field_main, bins, cbins, pdf, cg_wnoise, method)
                if projection == False: 
                    sigma_sq = None

                # OBTAIN P2D AND P3D           
                if 1:   
                    if lowpassfilter == 'no':
                        four_amps, powerspec2D, powerspec3D, infxyz, ks = powerspec.powerspec(cg_wnoise, dimension, target, many)    
                    if lowpassfilter == 'yes':
                        ks, powerspec2D, powerspec3D, infxyz = powerspec.powerspec_lowpass(cg_wnoise, the_k=ksnr[val], lowpass_filter=lowpassfilter)
                
                # OBTAIN SIGMA RHO
                if 1:
                    sigma_rho, sigma_rho2d, brunts_rsq = parsevals.parsevals_th(dimension, powerspec2D, powerspec3D, infxyz, ks, sigma_sq)
                    collect[target].append(brunts_rsq)   # what you are storing.


                
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

                rsqrt = np.append(rsqrt, the_rsallmed)
                if target == inf:
                    print('rsqrt pure ',the_rsallmed)
                #stop() #logmed ~ medlog, check [ ]

    # SAVE THE POWERSPECTRUM PLOTS 
    if 0:  
        plt.ylim(1e1,1e11)
        #plt.xlim(40,75)
        fig.savefig('pd2_snr_sn_wn_logy_%s_%d.png'%(box,proj_axis))

    # POSSIBLE FIT FUNCTIONS
    def cfunct(x, a, b, c, d):
        y = a * np.tanh((c*x)+d) + b
        return y

    def dfunct(x, a, b):
        y = 1/(1 + a*np.exp(-x+b)) 
        return y

    # for y-proj.
    rhalf_ctrue = 0.3035  #0.3035309747962126
    rhalf_ltrue = 0.3697  #0.3697364679188141
    rhalf_rtrue = 0.3385  #0.33851156031827845
    rhalf_btrue = 0.3274  #0.3273638724473809
    rhalf_ftrue = 0.3281  #0.3280843502158977
    # for x-proj.
    rhalf_ctrue_x =  0.34021165117508156
    # for z-proj.
    rhalf_ctrue_z = 0.36374026071471616
    if 0:
        # r^{1/2}true / r^{1/2}full vs SNR - raw
        ytrue_cbox = [0.30200493, 0.3006561,  0.30287845, 0.30259622, 0.30320516, 0.30334879, 0.30352811]
        ytrue_lbox = [0.36824365, 0.3592389,  0.36746576, 0.36903979, 0.37024448, 0.37046739, 0.3696961] 
        ytrue_rbox = [0.34350638, 0.3291993,  0.33969969, 0.33873154, 0.33918174, 0.33809821, 0.3384887]

        ydata = rsqrt[:-1]
        ydatalog = np.log10(ydata)
        xdata = target_snr[:-1]
        xdatalog = np.log10(xdata)
        if 0:  #tanh
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
        plt.ylabel(r'$R^{1/2}_{full}$') #plt.ylabel(r'$R^{1/2}_{true}/R^{1/2}_{full}$')
        plt.xlabel('SNR')
        #ax.set(ylim=[2e-1,1.1e0])  #ax.set(xlim=[1e-3,2e2], ylim=[1e-1,4e-1])
        plt.savefig('rsqrt_full_vs_SNR_%s_%d_%s'%(box,proj_axis,void_size))  #rsqrt_full, rsqrt_true, rhalf_ratio
        plt.clf()



    # to learn which are the crossing k values - then these get plotted at ''CALLING BOXPLOTS'' 
    print('okay')
    if 1:
        kcollect = defaultdict(list)  #CAREFUL?
        kmatch1 = np.linspace(1,16,16)  #target 0.01, 0.1
        kmatch2 = np.linspace(10,25,16)  #target 1
        kmatch3 = np.linspace(20,35,16)  #target:5,3
        kmatch4 = np.linspace(34,54,11)  #target 10
        kmatch5 = np.linspace(50,60,11)   #target 100
        iterations = 1000
        target_snr = target_snr[:-1]
        kmatch = [kmatch1, kmatch1, kmatch2, kmatch3, kmatch3, kmatch4, kmatch5]
        for i,k in enumerate(kmatch):
            kmatching = kmatch[i]
            target = target_snr[i]
            for val,ksnr in enumerate(kmatching):
                for iterate in range(iterations):  
                    cg_wnoise, noiseonly = noiseimage.add_noise(cg_first, target)
                    sigma_sq = calculate.Calculate(field_main, the_region, cg_wnoise, cg_second) 
                    ks, powerspec2D, powerspec3D, infxyz = powerspec.powerspec_lowpass(cg_wnoise, the_k=ksnr, lowpass_filter='yes')
                    sigma_rho, sigma_rho2d, brunts_rsq = parsevals.parsevals_th(dimension, powerspec2D, powerspec3D, infxyz, ks, sigma_sq)
                    kcollect[ksnr].append(brunts_rsq)   

            # CALLING BOXPLOTS; README 4 (an example). README 5 (plotting, dealing with BIG numbers in axis).
            # TRUE values - l_c_r: 0.269, 0.308, 0.265; b_c_f: 0.288, 0.308(yay), 0.243 
            if 1:  
                labels, data = kcollect.keys(), kcollect.values()  #PICK THE RIGHT DICTIONARY
                plt.boxplot(data, showfliers=True, showmeans=False)

                the_mean = rhalf_ctrue_x  #rsqrt at inf = pure signal - automate at least for the three projections: HERE!! 
                plt.plot([1,len(kmatch[i])],[the_mean,the_mean],'--',c='k',alpha=0.7)
                plt.xlabel('k')
                plt.ylabel(r'$R^{1/2}$')
                plt.xticks(range(1, len(labels) + 1), labels)
                #plt.yscale('log')
                #plt.ylim(0.10,0.39)
                plt.savefig('best_kvalue_%s_00%d_%d_%s_wpure'%(box,target,proj_axis,void_size))
                print(box,target)
                plt.clf()


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




# purpose 
'''
i. slice of positions of the cubes 
ii. plot 2d power spectrums of noise only and signal + noise

# for any box(position, proj, void) for each SNR 
a. get R^{1/2}_full values, where final is R^{1/2}_pure - call SNR function(k values turned off)

b. find the k values: AUTOMATE [ ] 
    1. based on ii, for a possible range of k values for a specific SNR target, boxplot R^{1/2} vs that k range; iterations = 1000
    2. pick the k value that matches R^{1/2} with R^{1/2}_inf most closely BY HAND; store 

c. get R^{1/2}_true values - call SNR function(k values turned on)

# figures for any box(position, proj, void) or boxes combined 'superplot'
# plot or boxplot ...
1. R^{1/2}_full vs SNR
2. R^{1/2}_true vs SNR with R^{1/2}_pure

3a. R^{1/2}_true/R^{1/2}_full vs SNR
3b. Fit 3a. - Make a table with all the parameters

4?. k values vs SNR
'''


