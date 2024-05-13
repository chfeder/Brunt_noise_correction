
from startlaser import * 
from collections import defaultdict
# for slices we need this encantation
from matplotlib import rcParams
rcParams['text.latex.preamble'] = r'\newcommand{\mathdefault}[1][]{}'

void_size = '50mu'
the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_%s_bubbles/NIF_hdf5_plt_cnt_0375'%(void_size))

# SETUP:
field_main = 'density'
field_two = 'cell_volume'
field_three = None 

proj_axis = 0  #or simply axis 
polytropic_index = 1 # +/-0.1, for 50um voids
the_level = 0

# CONSTRUCTING THE VOLUMES, pick x and z:
if 0:
    boxes = ['sx_box', 'sz_box']  
if 1:  #debug 
    boxes = ['sx_box']

# M-sized boxes.
delta_m = 0.02  #(1.2/3)/2 
edge_mtop = 0.3675 
edge_mbot = 0.3275
viz_top = 0.4
viz_bot = 0.0
# variable volumes
delta = 0.01


# FOR PLOTS 
if 1:  # MUST CLOSE BEFORE SECOND ITERATION
    fig,ax = plt.subplots(1,1)
    fig.subplots_adjust(wspace=0, hspace=0)
# POSSIBLE FITTING
def afunct(x, a, b):
    y = a * x + b
    return y

collect = defaultdict(list)
collect_3D = defaultdict(list)
for num, boxing in enumerate(boxes):  
    print('one iteration!')
    for pull in range(9): 
        pull = pull + 1
        if boxing == 'sx_box':
            the_ledge = the_ds.arr([-delta*pull, edge_mbot, -delta_m],'code_length') 
            the_redge = the_ds.arr([delta*pull, edge_mtop, delta_m],'code_length')  
            #the_ledge = the_ds.arr([-delta*pull, viz_bot, -delta_m],'code_length')  
            #he_redge = the_ds.arr([delta*pull, viz_top, delta_m],'code_length')  
        if boxing == 'sz_box':
            the_ledge = the_ds.arr([-delta_m, edge_mbot, -delta*pull],'code_length')
            the_redge = the_ds.arr([delta_m, edge_mtop, delta*pull],'code_length')
            #the_ledge = the_ds.arr([-delta_m, viz_bot, -delta*pull],'code_length')
            #the_redge = the_ds.arr([delta_m, viz_top, delta*pull],'code_length')

        the_center =0.5*the_ledge + 0.5*the_redge
        mindx = the_ds.index.get_smallest_dx()
        size = ((the_redge - the_ledge)/mindx).astype('int') 
        
        the_region = the_ds.region(the_center, the_ledge, the_redge) 
        section = extract_s1.extract_section_1(the_ds,the_level,the_ledge,size)
        
        # OBTAIN THE CUBE: 2D FIRST
        projection = True 
        dimension = 1
        cg_first,cg_second,cg_third = project.Project(section,proj_axis,projection,field_main,field_two,field_three)

        # OBTAIN SIGMA N^2
        print('getting information')
        sigma_sq, _ = calculate.Calculate(field_main, cg_first, cg_second) 
        if 0:
            collect[boxing].append(sigma_sq)
        if projection == False: 
            sigma_sq = None

        # OBTAIN P2D (AND P3D)           
        four_amps, powerspec2D, powerspec3D, infxyz, ks = powerspec.powerspec(cg_first, dimension)    
 
        # OBTAIN SIGMA RHO
        sigma_rho, sigma_rho2d, brunts_rsq, _ = parsevals.parsevals_th(dimension, powerspec2D, powerspec3D, infxyz, ks, sigma_sq)
        if 0:
            collect[boxing].append(brunts_rsq)   
        if 0: 
            collect[boxing].append(sigma_rho2d)   

        # BTW...
        if 0:  #HOW DOES IT LOOK 
            print('how does it look; slices')
            if pull == 1:
                a_slice = yt.SlicePlot(the_ds, proj_axis, ('gas', 'density'))#, data_source=the_region)
                outname = 'density_slice0375_%s_%d_%s.png'%(box, pull)
                a_slice.save(outname)
                stop()
        if 0:  # CHECK OUT THE POWER SPECTRUMS
            print('power spectrums')
            if pull == 2:  #delete when not plotting
                if dimension == 1:
                    cfp.plot(powerspec2D, ks, color='blue') 
                    cfp.plot(powerspec3D, ks, color='green', ylog=True,ylim=[10e1,10e10], save='powerspec3d2d_2d_%s_%s_%s'%(pull,box,void_size))
                if dimension == 2:
                    cfp.plot(powerspec3D, ks, color='orange', ylog=True,ylim=[10e1,10e10], save='powerspec3d_%s_%s_%s'%(pull,box,void_size))


        # OBTAIN THE CUBE: 3D NEXT
        print('getting the calculated version')
        projection = False 
        dimension = 2
        cg_first,cg_second,cg_third = project.Project(section,proj_axis,projection,field_main,field_two,field_three)
        _, sigma_sqrt = calculate.Calculate(field_main, cg_first, cg_second) 
        collect_3D[boxing].append(sigma_sqrt)   

        if 1:  #CROSS-CHECK WITH CF PIPELINE 
            dens_cube = cg_first  
            dlos = 1/dens_cube.shape[0]  
            coldens = turbulence_driving_analysis.get_column_density_map(dens_cube, dlos, 0) 
            coldens = coldens/coldens.mean() 
            power_spectrum_map, k = turbulence_driving_analysis.get_power_spectrum(coldens, zero_pad=False, weights=None, mirror=False, normalise=False, show='Fourier image')

            sigma_2D = np.std(coldens)
            Brunt_sqrt_R, power3D = turbulence_driving_analysis.get_Brunt_sqrt_R(power_spectrum_map, k)
            sigma_3D_from_Brunt = sigma_2D / Brunt_sqrt_R
            print('ready to collect cf pipe version')
            collect[boxing].append(sigma_rho2d)   

            sigma_3D_from_power = np.sqrt(power3D)
            sigma_3D_from_std = np.std(dens_cube)



    # PLOTS & FITS:  #LETS AUTOMATE THIS... 
    print('got the info, plotting time')
    the_d = np.linspace(0.02, 0.18, pull, endpoint=True)  
    the_l = delta_m * 2
    the_x = the_d/the_l # aspect ratio 
    the_xlog = np.log10(the_x)

    the_prey = [l for l in collect.values()]
    emptylist = []
    the_y = np.append(emptylist,the_prey[num])

    the_ylog = np.log10(the_y) 
    y_logerr_arr = np.array([0.046]*len(the_ylog))  #justify this number ..again 
    y_err_up = 10**(the_ylog + y_logerr_arr) - the_y 
    y_err_low = the_y - 10**(the_ylog - y_logerr_arr) 
    y_err_arr = [y_err_low,y_err_up] 

    the_preothery = [l for l in collect_3D.values()]
    emptylist_3D = []
    the_othery = np.append(emptylist_3D,the_preothery[num])

    # FITTING A CURVE
    print('will try to fit a curve now')
    if 1:
        rets = cfp.fit(afunct, the_xlog, the_ylog, yerr=y_logerr_arr)       
        the_xrecipe = np.linspace(the_xlog.min(),the_xlog.max(),num=500)  #try with np.logspace(), num=more numbers for other than a line, was num=len(the_xlog)
        the_xten =10**the_xrecipe
        the_yrecipe = afunct(the_xrecipe, *rets.popt)  #was the_yrecipe = (rets.popt[0]*the_xrecipe + rets.popt[1])
        the_yten = 10**the_yrecipe 
        if 1:
            cfp.plot(x=the_xten, y=the_yten, color='green', linestyle='dashed')
        print('param cf vals ',rets.popt)

    print('will try to fit a curve now')
    # PLOT THE CURVE; using cpack 
    if 1:
        boolean = True  #True: log axes, False: no log axes 
        if 1:
            cfp.plot(x=the_x, y=the_othery, type='scatter', color='orange', xlog=boolean, ylog=boolean)
            cfp.plot(x=the_x, y=the_y, yerr=y_err_arr, marker='o', linestyle='none', xlabel=r'$los/l$', ylabel=r'$\sigma_{\rho/\rho_o, BruntCF}$', \
                 xlog=boolean, ylog=boolean, save='sigma3d_Brunt_%s_aspect_logaxes_nolims_cfpipe_%s.png'%(boxing,void_size))
                 #ylabel=r'$\sigma_{\rho/\rho_o}$', r'$R^{1/2}$', r'$\sigma_{N/N_o}$', r'$\sigma_{Brunts}/\sigma_{3D}$'
    print('saving, saved')




 # DELETE??
if 0:
    if box == 'sx_box':
        if void_size == '50um':
            the_tsigmaR = [0.31001771892983476, 0.3078297687398883, 0.30420254190841384, 0.2995425654416588, 0.29491911865113646, 0.2943182444076409,
                           0.2937636729283891, 0.295444251435053, 0.3086992336937396]  
            the_tsigmaN = [0.16695969971020455, 0.1262302783976615, 0.10756128105834616, 0.10012803155043495, 0.09017580278138583, 0.09051612269382978,
                           0.08342417252845549,  0.08255569421346, 0.08220358374612546] 
        if void_size == '100um':
            the_tsigmaR = [0.3170498136892316, 0.34611830556097223, 0.3841385107122287, 0.385818914457437, 0.4084524344623558, 0.416317816355842,
                           0.41795659400450597, 0.4226154700364052, 0.4260310378576606]  
            the_tsigmaN = [0.15293937128752827, 0.11951510496293904, 0.10454399714497241, 0.09526700669557236, 0.11724506206874945, 0.13016455542450991,
                           0.12016730026984777, 0.11236920770172563, 0.09983649177817172] 
        if void_size == '12p5um':
            the_tsigmaR = [0.3342287396747619, 0.3178551701066717, 0.30462250105712285, 0.28594488288323616, 0.2698658383811071, 0.25892791745408256,
                           0.2557681991148804, 0.254177112629894, 0.2791405956929511] 
            the_tsigmaN = [0.2515456051989192, 0.2271953770726194, 0.20479582538332905,0.17250967118676916, 0.149931490568065, 0.13344083221736086,
                           0.12263421880655324, 0.11603307983897987, 0.1147919350666477]  
    if box == 'sz_box':
        if void_size == '50um':
            the_tsigmaR = [0.29975269098689805, 0.3078297687398883, 0.318353770910557, 0.3112775959416291, 0.30599288334855623, 0.3004886389624372,
                           0.29632498932085727, 0.2928859556167479, 0.2961103448111368]  
            the_tsigmaN = [0.1785115610856251,0.13833599403292984, 0.11492547915476435, 0.09961139448878231, 0.09145275071159353, 0.08284709980295706,
                           0.07515207490697227, 0.06966944077062515, 0.06643609685806992]
        if void_size == '100um':
            the_tsigmaR = [0.3576944256026828, 0.34611830556097223, 0.3750201357528904, 0.41319654960414226, 0.4349071434450878, 0.43806737183309685, 
                           0.4309897157698949, 0.4248607377476654,0.4234091800654715]  
            the_tsigmaN = [0.18859509606876412, 0.11697689397438824,0.10417725508647176,0.10754343915085687, 0.09995770608117312,0.09695294335454548,
                           0.0831420414386911, 0.07673816453855432, 0.0708454022364373]
        if void_size == '12p5um':
            the_tsigmaR = [0.32979540128257756,0.3178551701066717, 0.2976554884518944, 0.2805663944818849, 0.2665967792942205, 0.2553934628391015,
                           0.2500127693861195, 0.24635852020007531, 0.2646292996131689] 
            the_tsigmaN = [0.24406163181652032, 0.2240330567092556, 0.1976677896205206, 0.17475141863190619, 0.1519246484688996, 0.13493507695323137,
                           0.1265055393613463, 0.12023505563587977, 0.11868052152719089]


    if 0:  #lets leave on the side for now
        # R SQRT
        the_tRsqrt = np.array(the_tsigmaN)/np.array(the_tsigmaR)       
        the_othery = the_tsigmaR  
        the_simgaratios = the_yvalues/the_tsigmaR 