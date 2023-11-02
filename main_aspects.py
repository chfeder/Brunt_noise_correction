
from startlaser import * 
from collections import defaultdict

# for slices we need this encantation
from matplotlib import rcParams
rcParams['text.latex.preamble'] = r'\newcommand{\mathdefault}[1][]{}'

if 1: #12.5mu bubbles
    the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_12p5mu_bubbles/NIF_hdf5_plt_cnt_0375')
if 0: #50mu bubbles
    the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_50mu_bubbles/NIF_hdf5_plt_cnt_0375')
if 0: #100mu bubbles
    the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_100mu_bubbles/NIF_hdf5_plt_cnt_0375')

'''
to obtain SIGMA_RHO:
'''
# SETUP:
field_main = 'density'
field_two = 'cell_volume'
field_three = None 

projection = True 
dimension = 1
noise = False
noiseplot_only = False

proj_axis = 2  #or simply axis 
polytropic_index = 1 # +/-0.1, for 50um voids
the_level = 0

# CONSTRUCTING THE VOLUMES, pick x or z:
box = 'sz_box'  

# M-sized boxes.
delta_m = 0.02  #(1.2/3)/2 
edge_mtop = 0.3675 
edge_mbot = 0.3275
viz_top = 0.4
viz_bot = 0.0

# variable volumes
delta = 0.01

# SIGNAL TO NOISE RATIOs: 
target_snr = [inf]

iterations = 1
noiseat = 'none'  # when to add noise, at 'image' or at 'dft - no noise really...'
lowpassfilter = 'no'  # yes:low k values, no:all k values

method = 0  # method to obtain sigma N - 0: via PDF, 1: via np.std(), 2: rms 
imaging = 'no'  # make images?
width = .12  # for images


collect = defaultdict(list)
# FOR PLOTS 
if 1:
    fig,ax = plt.subplots(1,1)
    fig.subplots_adjust(wspace=0, hspace=0)

for pull in range(9): 
    pull = pull + 1
    if box == 'sx_box':
        the_ledge = the_ds.arr([-delta*pull, edge_mbot, -delta_m],'code_length') 
        the_redge = the_ds.arr([delta*pull, edge_mtop, delta_m],'code_length')  
        #the_ledge = the_ds.arr([-delta*pull, viz_bot, -delta_m],'code_length')  
        #he_redge = the_ds.arr([delta*pull, viz_top, delta_m],'code_length')  
    if box == 'sz_box':
        the_ledge = the_ds.arr([-delta_m, edge_mbot, -delta*pull],'code_length')
        the_redge = the_ds.arr([delta_m, edge_mtop, delta*pull],'code_length')
        #the_ledge = the_ds.arr([-delta_m, viz_bot, -delta*pull],'code_length')
        #the_redge = the_ds.arr([delta_m, viz_top, delta*pull],'code_length')

    the_center =0.5*the_ledge + 0.5*the_redge
    mindx = the_ds.index.get_smallest_dx()
    size = ((the_redge - the_ledge)/mindx).astype('int') 
    the_region = the_ds.region(the_center, the_ledge, the_redge) 

    section = extract_s1.extract_section_1(the_ds,the_level,the_ledge,size)
    section_blur = motion_blur.Motion_Blur()  #tbd; simply returns
    cg_first,cg_second,cg_third = project.Project(section,proj_axis,projection,field_main,field_two,field_three)

    if 0:  
        if pull == 1:
            a_slice = yt.SlicePlot(the_ds, proj_axis, ('gas', 'density'))#, data_source=the_region)
            outname = 'density_slice0375_%s_%d_%s.png'%(box, proj_axis, pull)
            a_slice.save(outname)
            stop()

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
                #print(pull)
                #print(sigma_rho)


# PLOTS & FITS: 
if 1:
    if box == 'sx_box':
        #the_truey = [0.31001771892983476, 0.3078297687398883, 0.30420254190841384, 0.2995425654416588, 0.29491911865113646, 0.2943182444076409,
        #             0.2937636729283891, 0.295444251435053, 0.3086992336937396] #50mu
        #the_truey = [0.3170498136892316, 0.34611830556097223, 0.3841385107122287, 0.385818914457437, 0.4084524344623558, 0.416317816355842,
        #             0.41795659400450597, 0.4226154700364052, 0.4260310378576606]  #100mu
        the_truey = [0.3342287396747619, 0.3178551701066717, 0.30462250105712285, 0.28594488288323616, 0.2698658383811071, 0.25892791745408256,
                     0.2557681991148804, 0.254177112629894, 0.2791405956929511]  #12p5mu
    if box == 'sz_box':
        #the_truey = [0.29975269098689805, 0.3078297687398883, 0.318353770910557, 0.3112775959416291, 0.30599288334855623, 0.3004886389624372,
        #             0.29632498932085727, 0.2928859556167479, 0.2961103448111368]  #50um
        #the_truey = [0.3576944256026828, 0.34611830556097223, 0.3750201357528904, 0.41319654960414226, 0.4349071434450878, 0.43806737183309685, 
        #             0.4309897157698949, 0.4248607377476654,0.4234091800654715]  #100mu
        the_truey = [0.32979540128257756,0.3178551701066717, 0.2976554884518944, 0.2805663944818849, 0.2665967792942205, 0.2553934628391015,
                     0.2500127693861195, 0.24635852020007531, 0.2646292996131689]  #12p5mu
           
    the_d = np.linspace(0.02, 0.18, pull, endpoint=True)  
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
    if 1:
        rets = cfp.fit(afunct, the_xlog, the_ylog, yerr=y_logerr_arr)
        
        the_xrecipe = np.linspace(the_xlog.min(),the_xlog.max(),num=500)  #try with np.logspace(), num=more numbers for other than a line, was num=len(the_xlog)
        the_xten =10**the_xrecipe
        the_yrecipe = afunct(the_xrecipe, *rets.popt)  #was the_yrecipe = (rets.popt[0]*the_xrecipe + rets.popt[1])
        the_yten = 10**the_yrecipe
        
        cfp.plot(x=the_xten, y=the_yten, color='green', linestyle='dashed')
        print('param cf vals ',rets.popt)
    # PLOT THE CURVE; using cpack 
    if 1:
        cfp.plot(x=the_x, y=the_truey, type='scatter', color='orange')
        cfp.plot(x=the_x, y=the_y, yerr=y_err_arr, marker='o', linestyle='none', xlabel='aspect ratio', \
                 ylabel=r'$\sigma_{\rho/\rho_o}$', xlog=True, ylog=True, save='sigma_rho_%s_aspect_logaxes_nolims_fit_errbars_wcfp_12p5um.png'%(box))
#ylim=[1.8e-1,5.5e-1]
