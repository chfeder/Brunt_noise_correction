
from startlaser import * 
from collections import defaultdict
# for slices we need this encantation
from matplotlib import rcParams
rcParams['text.latex.preamble'] = r'\newcommand{\mathdefault}[1][]{}'
import brunt_tools_cp as bt
reload(bt)
import fourier_tools_py3.fourier_filter as Filter

void_size = '50mu'
the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_%s_bubbles/NIF_hdf5_plt_cnt_0375'%(void_size))

'''
to obtain SIGMA_RHO:
'''
# SETUP:
field_main = 'density'
field_two = 'cell_volume'
field_three = None 

polytropic_index = 1 # +/-0.1, for 50um voids
the_level = 0

# CONSTRUCTING THE VOLUMES, pick x or z:
box = 'sx_box'  
proj_axis = 0  #proj_axis needs to match  

# M-sized boxes.
delta_m = 0.02  #(1.2/3)/2 
edge_mtop = 0.3675 
edge_mbot = 0.3275
viz_top = 0.4
viz_bot = 0.0

# variable volumes
delta = 0.01


the_y_2d = []
the_y_3d2d = []
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
    section = extract_s1.extract_section_1(the_ds,the_level,the_ledge,size)

    the_region = 0

    if pull == 2:
        if 1:  #get your 3D cube.
            projection = False 
            dimension = 2
            cg_first_3d,_,_ = project.Project(section, proj_axis, projection, field_main) 
            #_, p2d, p3d, dk3, ks = powerspec.powerspec(cg_first_3d, dimension)    
            #sigma_rho, sigma_rho2d, brunts_rsq = parsevals.parsevals_th(dimension, p2d, p3d, dk3, ks)
            #print("sig rho, sig rho2d ", sigma_rho, sigma_rho2d)
            #cfp.plot(p3d, ks, color='orange', ylog=True,ylim=[10e1,10e10], save='powerspec3d_3d2d_%s_%s_%s'%(pull,box,void_size))

        if 0:  #sigma brunt using dcc code 
            ftool = bt.fft_tool(cg_first_3d, cg_second)
            ftool.do3()
            ftool.do2(projax=proj_axis)
            sig_B, sigx2d, rsqrt = bt.sigmas_back(ftool, dimension, method='norm')  

            print("sigma_2D dcc_etal = ", sigx2d)
            print("Brunt R^(1/2) dcc_etal = ", rsqrt)
            print("sigma_3D_from_Brunt dcc_etal = ", sig_B)
            print("dimension", dimension)

        if 0:  #sigma brunt using ljv-dcc code
            projection = True 
            dimension = 1
            cg_first,cg_second,_ = project.Project(section, proj_axis, projection, field_main, field_two) 
            sigma_sq = calculate.Calculate(field_main, the_region, cg_first, cg_second) 
            _, p2d, p32d, dk32, ks = powerspec.powerspec(cg_first, dimension)    
            sigma_rho, sigma_rho2d, brunts_rsq, sigmaN = parsevals.parsevals_th(dimension, p2d, p32d, dk32, ks, sigma_sq)
            '''
            print("sigma_2D_ljv_dcc = ", np.sqrt(sigmaN))
            print("Brunt R^(1/2) ljv_dcc = ", brunts_rsq)
            print("sigma_3D_from_Brunt ljv_dcc = ", sigma_rho)
            print("        ")
            cfp.plot(p2d, ks, color='blue') #,ylog=True,ylim=[10e1,10e11], save='powerspec3d2d_%s_%s'%(pull,box))
            cfp.plot(p32d, ks, color='green') #,ylog=True,ylim=[10e1,10e11], save='powerspec3d2d_%s_%s'%(pull,box)) 
            '''

        if 0:  #sigma brunt using cf_ig code 
            dens_cube = cg_first_3d  
            dlos = 1/dens_cube.shape[0]  #check why 0... 
            #dens = np.concatenate([dens_cube]*1, 1)
            # compute column density
            coldens = turbulence_driving_analysis.get_column_density_map(dens_cube, dlos, 0) 
            coldens = coldens/coldens.mean() 
            power_spectrum_map, k = turbulence_driving_analysis.get_power_spectrum(coldens, zero_pad=False, weights=None, mirror=False, normalise=False, show='Fourier image')

            sigma_2D = np.std(coldens)
            Brunt_sqrt_R, power3 = turbulence_driving_analysis.get_Brunt_sqrt_R(power_spectrum_map, k)
            sigma_3D_from_Brunt = sigma_2D / Brunt_sqrt_R
            sigma_3D_from_power = np.sqrt(power3)
            sigma_3D_from_std = np.std(dens_cube)
            the_y_2d = np.append(the_y_2d, sigma_3D_from_Brunt)
            the_y_3d2d = np.append(the_y_3d2d, sigma_3D_from_power)

            #print("sigma_2D cf_ig = ", sigma_2D)
            #print("Brunt R^(1/2) cf_ig = ", Brunt_sqrt_R)
            #print("sigma_3D_from_Brunt cf_ig= ", sigma_3D_from_Brunt)
            #print("        ")

        # POWER SPECTRUM PLOTS with cg_ig power_spectrum_map
        if 0:
            # using ljv_dcc code. 
            power_spectrum_flat = power_spectrum_map.flatten()
            kflat = k.flatten()
            kbins = np.arange(0, kflat.max()+1, 1) 
            kvals = 0.5 * (kbins[1:] + kbins[:-1])
            Abins, _, _ = stats.binned_statistic(kflat, power_spectrum_flat, \
                          statistic = "sum", bins = kbins)  
            cfp.plot(Abins, kvals, color='blue', ylog=True, xlog=False, save='ps_2dB_cf_ljv%s_%s_%s'%(pull,box,void_size)) 
        if 1:
            # using cfpack
            stop()
            rets = cfp.get_spectrum(cg_first_3d) #power_spectrum_map
            cfp.plot(Abins, kvals, color='blue', ylog=True, xlog=False, save='ps_2dB_cf_ljv%s_%s_%s'%(pull,box,void_size)) 

        if 0:
            # plot the power spectrum using dcc_etal code. 
            print('hello')
            power_spectrum_map/=power_spectrum_map.size  #WHY
            #stop()
            ff2 = Filter.FourierFilter(power_spectrum_map)
            power_1d2 = np.array([power_spectrum_map[ff2.get_shell(bin)].sum() for bin in range(ff2.nx)])
            k2d=ff2.get_shell_k()
            cfp.plot(power_1d2, k2d, color='blue', ylog=True, xlog=False, save='ps_2dB_cf_dcc%s_%s_%s'%(pull,box,void_size)) 



        # PREVIOUSLY BOTHCODES.PY
        if 0:
            proj_axis = 0
            projection = True 
            dimension = 1
            cg_first, cg_second,_ = project.Project(section, proj_axis, projection, field_main, field_two) 
            _, p2d, p32d, dk32, ks = powerspec.powerspec(cg_first, dimension)    

            projection = False 
            dimension = 2
            cg_first_3D, cg_second,_ = project.Project(section, proj_axis, projection, field_main, field_two) 
            _, p2dd, p32dd, dk32, kss = powerspec.powerspec(cg_first_3D, dimension)    

            cg_norm_2d = cg_first/cg_first.mean() - 1 
            cg_norm_3d = cg_first_3D/cg_first_3D.mean() - 1 
            ftool = bt.fft_tool(None)
            ftool.rho2 = cg_norm_2d
            ftool.do2(projax=proj_axis)


            # PLOTTINGS
            fig,ax=plt.subplots(1,1)
            # script1: ljv_dcc, script2: dcc, script 3: cf
            ax.plot(ks, p2d/p2d[3], c='b', label='2D_script1')
            ax.plot(ftool.k2d/ftool.k2d[1], ftool.power_1d2/ftool.power_1d2[3], c='b', linestyle='dashed', label='2D_script2')

            ftool = bt.fft_tool(cg_norm_3d)
            ftool.do3()
            ax.plot(ks, p32d/p32d[3], c='g', label='3D2D_script1')
            ax.plot(kss, p32dd/p32dd[3], c='g', linestyle='dotted', label='3D2D_script1') #what is the difference here
            ax.plot(ftool.k3d/ftool.k3d[1], ftool.power_1d3/ftool.power_1d3[3], c='g', linestyle='dashed', label='2D_script2')

            ax.set(yscale='log')
            ax.set(ylim = [10e-9,10e2])
            fig.savefig('power23_twoscripts')
            
            print('saved plots')        

# before PREVIOUSLY BOTHCODES.PY
'''
if 0:
    the_d = np.linspace(0.02, 0.18, pull, endpoint=True)  
    the_l = 0.02 * 2
    the_x = the_d/the_l  #this is the aspect ratio 

    boolean = False
    fig,ax=plt.subplots(1,1)
    ax.scatter(the_x, the_y_3d2d, color='orange')
    ax.scatter(the_x, the_y_2d, color='blue')
    fig.savefig('sigma_brunt_cf_ig_w3d2d_%s_%s'%(box,void_size))
'''


        




