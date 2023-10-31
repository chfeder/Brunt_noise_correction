
from startlaser import * 
from collections import defaultdict
the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_50mu_bubbles/NIF_hdf5_plt_cnt_0375')

'''
to obtain SIGMA_RHO:
'''
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

# CONSTRUCTING THE VOLUMES, pick x or z:
box = 'sx_box'  

# M-sized boxes.
delta_m = 0.02  #(1.2/3)/2 
edge_mtop = 0.3675 
edge_mbot = 0.3275

# variable volumes
delta = 0.01

# SIGNAL TO NOISE RATIOs: 
target_snr = [inf]

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
                four_amps, powerspec2D, powerspec3D, infxyz, ks = powerspec.powerspec(cg_wnoise, dimension, box, target, many)    
            if 0:  # if adding noise to DFT: DELETE!
                noisy_powerspec2D, noisy_powerspec3D, noise_sign, noise_signs = signaltonoise.signaltonoise(cg_first, dimension, target, four_amps, powerspec2D, powerspec3D) 
            
            # OBTAIN SIGMA RHO
            if 1:
                sigma_rho, brunts_rsq = parsevals.parsevals_th(dimension, powerspec2D, powerspec3D, infxyz, ks, sigma_sq)
                #collect[target].append(brunts_rsq)   #plotting brunts_rsq  
                collect[target].append(sigma_rho)    #plotting sigma_rho 
                print(pull)
                print(sigma_rho)


# PLOTS & FITS: 
if 1:
    if box == 'sx_box':
        the_truey = [0.31001771892983476, 0.3078297687398883, 0.30420254190841384, 0.2995425654416588, 0.29491911865113646, 0.2943182444076409]
    if box == 'sz_box':
        the_truey = [0.29975269098689805, 0.3078297687398883, 0.318353770910557, 0.3112775959416291, 0.30599288334855623, 0.3004886389624372]

    the_d = np.linspace(0.02, 0.12, pull, endpoint=True)  
    the_l = delta_m * 2
    the_x = the_d/the_l # aspect ratio 
    the_xlog = np.log10(the_x)

    the_prey = [l for l in collect.values()]
    emptylist = []
    the_y = np.append(emptylist,the_prey[0])
    the_ylog = np.log10(the_y) 

    y_logerr_arr = np.array([0.046]*len(the_ylog)) 
    y_err_up = 10**(the_ylog + y_logerr_arr) - the_y 
    y_err_low = the_y - 10**(the_ylog - y_logerr_arr) 
    y_err_arr = [y_err_low,y_err_up] 


    # POSSIBLE FITTING
    def afunct(x, a, b):
        y = a * x + b
        return y
    # FITTING A CURVE
    if 0:
        rets = cfp.fit(afunct, the_xlog, the_ylog, yerr=y_logerr_arr)
        
        the_xrecipe = np.linspace(the_xlog.min(),the_xlog.max(),num=500)  #try with np.logspace(), num=more numbers for other than a line, was num=len(the_xlog)
        the_xten =10**the_xrecipe
        the_yrecipe = afunct(the_xrecipe, *rets.popt)  #was the_yrecipe = (rets.popt[0]*the_xrecipe + rets.popt[1])
        the_yten = 10**the_yrecipe
        
        cfp.plot(x=the_xten, y=the_yten, color='green', linestyle='dashed')
        print('param cf vals ',rets.popt)
    # PLOT THE CURVE; using cpack 
    if 0:
        cfp.plot(x=the_x, y=the_truey, type='scatter', color='orange')
        cfp.plot(x=the_x, y=the_y, ylim=[1.8e-1,5.5e-1], yerr=y_err_arr, marker='o', linestyle='none', xlabel='aspect ratio', \
                 ylabel=r'$\sigma_{\rho/\rho_o}$', xlog=True, ylog=True, save='sigma_rho_%s_aspect_logaxes_fit_errbars_wcfp.png'%(box))

