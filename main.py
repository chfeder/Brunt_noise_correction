
from startlaser import * 
from collections import defaultdict
# README: 0 (info on data sets)
the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_50mu_bubbles/NIF_hdf5_plt_cnt_0375')

'''
to obtain SIGMA_RHO:
if field_main = any velocity, then field_two = 'density', field_three = pressure
README: 1 (field considered)
'''

# SETUP:
field_main = 'density'
field_two = 'cell_volume'
field_three = None 

projection = True
dimension = 1
noise = False
noiseplot_only = False
noisy_powerspec2D = 0
noisy_powerspec3D = 0

proj_axis = 2  #or simply axis 
polytropic_index = 1 # +/-0.1, for 50um voids
the_level = 0

# which boxes
#boxes = ['l_box', 'c_box', 'r_box']
boxes = ['sz_box']

# "small" boxes
delta_a = 0.0175 
edge_atop = 0.365 
edge_abot = 0.33
# "big" boxes ...too big...
delta_b = 0.02667  #(1.6/3)/2 
edge_btop = 0.37417 
edge_bbot = 0.32083
# "medium" boxes
delta_m = 0.02  #(1.2/3)/2 
edge_mtop = 0.3675 
edge_mbot = 0.3275
# "stretchy box"
delta = 0.01

width = 2*delta_m  #for images

# README: 2 (in decibels)
target_snr = [inf]
#target_snr = [5, 10, 100]
#target_snr = [1, 3, 5, 10, 100, inf]
#target_snr = [0.001, 0.01, 0.1, 1, 3, 5, 10, 100, inf]
#target_colors = ['red','orange','y','cyan','magenta','blue','green','grey','black']
target_colors = ['blue','green','orange']
iterations = 1
method = 0  #method to obtain sigma N - 0: via PDF, 1: via np.std(), 2: rms 
imaging = 'no'  #make images?
noiseat = 'dft'  #when to add noise, at 'image' or at 'dft'


# THREE BOXES
# limit = np.zeros(4) #for README 5
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

    #THERE'S GOTTA BE A BETTER WAY TO DO THIS: SX & SZ BOXING.
    collect = defaultdict(list)
    for pull in range(6): 
        pull = pull + 1
        if box == 'sx_box':
            the_ledge = the_ds.arr([-delta*pull, edge_mbot, -delta_m],'code_length')
            the_redge = the_ds.arr([delta*pull, edge_mtop, delta_m],'code_length')
        if box == 'sz_box':
            the_ledge = the_ds.arr([-delta_m, edge_mbot, -delta*pull],'code_length')
            the_redge = the_ds.arr([delta_m, edge_mtop, delta*pull],'code_length')

    # indentation for the sx and sz boxes...
        the_center =0.5*the_ledge + 0.5*the_redge
        mindx = the_ds.index.get_smallest_dx()
        size = ((the_redge - the_ledge)/mindx).astype('int') 
        the_region = the_ds.region(the_center,the_ledge, the_redge)

        section = extract_s1.extract_section_1(the_ds,the_level,the_ledge,size)
        section_blur = motion_blur.Motion_Blur()  #tbd; simply returns
        cg_first,cg_second,cg_third = project.Project(section,proj_axis,projection,field_main,field_two,field_three)

        if 0:   
            images = make_images.Make_Images(the_ds,cg_first,proj_axis,box,the_center,width) 
        if 0:
            mom_map = momentmap.moment_map(cg_first)  #a brief exploration, maybe discard
        
        if 1:
            #collect = defaultdict(list)  #CAREFUL
            # p2Ds
            fig,ax = plt.subplots(1,1)
            fig.subplots_adjust(wspace=0, hspace=0)
            for val,target in enumerate(target_snr):
                for many in range(iterations):  

                    # ADD NOISE TO THE IMAGE
                    if 1:  
                        if noiseat == 'image':
                            cg_wnoise, cg_wnoisetwo, noiseonly = noiseimage.noiseimage(cg_first, cg_second, target, box, imaging)
                        if noiseat == 'dft':  #also works for nonoise
                            cg_wnoise = cg_first
                            cg_wnoisetwo = cg_second
                            noiseonly = 0
            
                    # OBTAIN SIGMA N^2
                    if 1:  #calculated
                        sigma_sq = calculate.Calculate(field_main, the_region, cg_wnoise, cg_wnoisetwo) 
                    if 0:  #from PDF
                        # the following could be split into two definitions/files PDF and FIT
                        bins, cbins, pdf = make_pdfNfit.Make_PDFnFIT(field_main, projection, proj_axis, cg_wnoise, cg_wnoisetwo,cg_third) 
                        sigma_sq = make_sigmasquared.Make_SigmaSquared(field_main, bins, cbins, pdf, cg_wnoisetwo, method)
                    if projection == False: 
                        sigma_sq = None

                    # OBTAIN P2D AND P3D           
                    if 1:   
                        four_amps, powerspec2D, noisepowerspec2D, powerspec3D, infxyz, ks, ksnoise = powerspec.PowerSpec(cg_wnoise, cg_wnoisetwo, noiseonly, noiseplot_only, dimension, box, target, many)
                        if 1:  # power spectrum plots
                            #plt.plot(ks, powerspec2D, c=target_colors[val], alpha=0.6)
                            if noiseplot_only == True:
                                plt.plot(ksnoise, noisepowerspec2D, c=target_colors[val], alpha=0.6)
                    if 0:  # if adding noise to DFT 
                        noisy_powerspec2D, noisy_powerspec3D, noise_sign, noise_signs = signaltonoise.signaltonoise(cg_first, dimension, target, four_amps, powerspec2D, powerspec3D) 
                        # README 3. (thoughts on making pdfs)
                    
                    # OBTAIN SIGMA RHO
                    if 1:
                        sigma_rho, brunts_rsq = parsevals.parsevals_th(dimension, powerspec2D, powerspec3D, infxyz, ks, noiseat, sigma_sq, noise, noisy_powerspec2D, noisy_powerspec3D)  
                        #collect[target].append(brunts_rsq)   #plotting brunts_rsq  
                        collect[target].append(sigma_rho)    #plotting sigma_rho 


    # PLOTS & FITS
    if 1:
        the_x = np.linspace(0.02, 0.12, pull, endpoint=True)
        the_xlog = np.log10(the_x)
        the_prey = [l for l in collect.values()]
        emptylist = []
        the_y = np.append(emptylist,the_prey[0])
        the_ylog = np.log10(the_y) 
        the_ylogsem = stats.sem(the_ylog)
        the_ylogstd = stats.sem(the_ylog)

        # FITTING A CURVE
        def afunct(x, a, b):
            #y = a*np.exp(b*x)
            #y = -a*np.log(x)+b
            #y = a*x + b
            #y = a * 1/np.sqrt(x)  
            y = a*x + b
            return y
        def bfunct(x, a, b):
            y = -0.5*x + a
            return y
        if 1:
            param, param_cov = curve_fit(afunct, the_xlog, the_ylog)
            plt.plot(the_xlog, afunct(the_xlog, *param), 'r--')
        if 1:
            param, param_cov = curve_fit(bfunct, the_xlog, the_ylog)
            plt.plot(the_xlog, bfunct(the_xlog, *param), 'b--')


        if box == 'sx_box':
            the_truey = [0.31001771892983476, 0.3078297687398883, 0.30420254190841384, 0.2995425654416588, 0.29491911865113646, 0.2943182444076409]
            the_trueylog = np.log10(the_truey) 
        if box == 'sz_box':
            the_truey = [0.29975269098689805, 0.3078297687398883, 0.318353770910557, 0.3112775959416291, 0.30599288334855623, 0.3004886389624372]
            the_trueylog = np.log10(the_truey) 

        ax.errorbar(the_xlog, the_ylog, yerr=the_ylogstd, fmt='o')
        plt.scatter(the_xlog, the_trueylog, color='orange')

        #ax.set_yscale('log')
        #ax.set_xscale('log')
        plt.xlabel('codelength')
        plt.ylabel(r'$\sigma_{\rho/\rho_o}$')
        plt.savefig('sigma_rho_%s_wtrue_wfitlogoversqrts_err_.png'%(box))
        #stop()

    # PLOT POWER SPECTRUMS  
    if 0:
        ax.set_yscale('log')
        plt.xlabel('$k$')
        plt.ylabel('$P(k)$')
        plt.ylim(1e-2,1e13)
        #plt.xlim(0.8,1e2)
        plt.tight_layout()
        #plt.savefig('brunts_P2D_%s_%.3f_%d_loglog.png'%(box_loc,snr,rep))
        #plt.savefig('brunts_P2D_%s_logy.png'%(box))
        plt.savefig('brunts_P2D_%s_logy_noiseonly.png'%(box))

    # CALLING BOXPLOTS; README 4 (an example). README 5 (plotting, dealing with BIG numbers in axis).
    # TRUE values - l_c_r: 0.269, 0.308, 0.265; b_c_f: 0.288, 0.308(yay), 0.243 
    if 0:  
        labels, data = collect.keys(), collect.values()
        plt.boxplot(data, showfliers=False, showmeans=True)
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
        plt.ylim(0.10,0.39)
        plt.savefig('boxplot_1000n_snr_%s_xpos_image_bruntsrsqVSsnr_001_fixed'%box)
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

