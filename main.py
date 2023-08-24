
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
noise = True
noisy_powerspec2D = 0
noisy_powerspec3D = 0
proj_axis = 1  #or simply axis 
polytropic_index = 1 # +/-0.1, for 50um voids
the_level = 0

# which boxes
boxes = ['b_box', 'c_box', 'f_box']
colors = ['green','blue','orange']

# X POSITION
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

width = 2*delta_m  #for images

# README: 2 (in decibels)
target_snr = [1,3,5,10,100, inf]
iterations = 1000
method = 0  #method to obtain sigma N - 0: via PDF, 1: via np.std(), 2: rms 
imaging = 'no'  #make images?


# THREE BOXES
# limit = np.zeros(4) #for README 5
for num, box in enumerate(boxes): 
    if box == 'c_box':
        the_ledge = the_ds.arr([-delta_m, edge_mbot, -delta_m],'code_length')
        the_redge = the_ds.arr([delta_m, edge_mtop, delta_m],'code_length')
    #if box == 'r_box':
    #    the_ledge = the_ds.arr([delta_m, edge_mbot, delta_m],'code_length')
    #    the_redge = the_ds.arr([delta_m+2*delta_m, edge_mtop, delta_m+2*delta_m],'code_length')
    if box == 'r_box':
        the_ledge = the_ds.arr([delta_m, edge_mbot, -delta_m],'code_length')
        the_redge = the_ds.arr([delta_m+2*delta_m, edge_mtop, delta_m],'code_length')
    if box == 'f_box':
        the_ledge = the_ds.arr([-delta_m, edge_mbot, delta_m],'code_length')
        the_redge = the_ds.arr([delta_m, edge_mtop, delta_m+2*delta_m],'code_length')
    #if box == 'l_box':
    #    the_ledge = the_ds.arr([-delta_m-2*delta_m, edge_mbot, -delta_m-2*delta_m],'code_length')
    #    the_redge = the_ds.arr([-delta_m, edge_mtop, -delta_m],'code_length')
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
        collect = defaultdict(list)
        for target in target_snr:
            for many in range(iterations):  

                # ADD NOISE TO THE IMAGE
                if 1:  
                    imagewnoise = noiseimage.noiseimage(cg_first, target, box, imaging)
                    cg_wnoise = imagewnoise
                
                # OBTAIN SIGMA N^2
                if 1:  #calculated
                    sigma_sq = calculate.Calculate(field_main, the_region, cg_wnoise, cg_second) 
                if 0:  #from PDF
                    # the following could be split into two definitions/files PDF and FIT
                    bins, cbins, pdf = make_pdfNfit.Make_PDFnFIT(field_main, projection, proj_axis, cg_wnoise, cg_second,cg_third) 
                    sigma_sq = make_sigmasquared.Make_SigmaSquared(field_main, bins, cbins, pdf, cg_wnoise, method)
                if projection == False: 
                    sigma_sq = None

                # OBTAIN P2D AND P3D           
                if 1:   
                    four_amps, powerspec2D, powerspec3D, infxyz, ks = powerspec.PowerSpec(cg_wnoise, cg_second, dimension)
                if 0:  #previously added noise to DFT 
                    noisy_powerspec2D, noisy_powerspec3D, noise_sign, noise_signs = signaltonoise.signaltonoise(cg_first, dimension, target, four_amps, powerspec2D, powerspec3D) 
                    # README 3. (thoughts on making pdfs)
                
                # OBTAIN SIGMA RHO
                if 1:
                    sigma_rho = parsevals.parsevals_th(dimension, powerspec2D, powerspec3D, infxyz, ks, sigma_sq, noise, noisy_powerspec2D, noisy_powerspec3D) 
                    collect[target].append(sigma_rho)     


        # CALLING BOXPLOTS; README 4 (an example). README 5 (plotting, dealing with BIG numbers in axis).
        if 1:  
            labels, data = collect.keys(), collect.values()
            plt.boxplot(data,showmeans=True)
            '''
            if box == 'l_box':
                plt.plot([1,6],[0.238,0.238],'--',c=[0.5]*4)
            if box == 'c_box':
                plt.plot([1,6],[0.308,0.308],'--',c=[0.5]*4)
            if box == 'r_box':
                plt.plot([1,6],[0.236,0.236],'--',c=[0.5]*4)
            '''
            plt.xlabel('SNR')
            plt.ylabel(r'$\sigma_{\rho/\rho_o} $')
            plt.xticks(range(1, len(labels) + 1), labels)
            #plt.ylim(2.75,3.75)
            plt.savefig('boxplot_1000n_snr_%s_wtrue_medN_Zpos'%box)
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

