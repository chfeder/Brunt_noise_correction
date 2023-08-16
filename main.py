
from startlaser import * 
from collections import defaultdict
'''
main.py: in progress 

THE DATA: the simulations are in /scratch/ek9/ccf100/nif
* 9 different simulations with a few snapshots:
* NIF_hdf5_plt_cnt_0*
* frames = [000,125,250,275,300,325,350,375,425]
'''

the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_50mu_bubbles/NIF_hdf5_plt_cnt_0375')

'''
to obtain SIGMA_RHO:
if field_main = any velocity, then field_two = 'density', field_three = pressure
derived field (an option): to take the mean out, aka k=0 mode
'''

def add_field(thing):
    def overdensity(field,data):
        rho_0 = data.get_field_parameter('rho_0')
        #rho_ratio = (data['density']/rho_0).v
        #return (rho_ratio-1)
        # OR SAME THING:
        return(((data['density'] - rho_0)/rho_0))
    thing.add_field(YT_overdensity,function=overdensity)#,units='dimensionless')
add_field(the_ds)


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

# target signal to noise suggested ratios: inf,100(20db),10(10db), 5(6.9897db), 3(4.7712db), 1(0db)
#target_snr_db = [-100, -90, -80, -70, -60, -50, -40, -30, -20, -10, 0, 4.7712,6.9897, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
#target_snr = [0.0000000001,0.000000001,0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1,1,3,5,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000,10000000000]
target_snr_db = [0, 4.7712,6.9897, 10, 20]
target_snr = [inf]  #[1,3,5,10,100, inf]

# which box? right, left, center, not a box
#cube = 'c_box'
boxes = ['l_box', 'c_box', 'r_box']
colors = ['green','blue','orange']
delta_i = 0.0175 
delta = 0.02667  #(1.6/3)/2 
width = 2*delta  #for images


# THREE BOXES 
limit = np.zeros(4)
#fig,ax = plt.subplots(1,1)
#fig.subplots_adjust(wspace=0, hspace=0)
for num, box in enumerate(boxes):  
    if boxes[num] == 'c_box':
        the_ledge = the_ds.arr([-delta,0.32083,-delta],'code_length')
        the_redge = the_ds.arr([delta,0.37417,delta],'code_length')
    if boxes[num] == 'r_box':
        the_ledge = the_ds.arr([delta,0.32083,delta],'code_length')
        the_redge = the_ds.arr([delta+2*delta,0.37417,delta+2*delta],'code_length')
    if boxes[num] == 'l_box':
        the_ledge = the_ds.arr([-delta-2*delta,0.32083,-delta-2*delta],'code_length')
        the_redge = the_ds.arr([-delta,0.37417,-delta],'code_length')
    '''
    if boxes[num] == 'notabox':
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
        box_loc = boxes[num]
        images = make_images.Make_Images(the_ds,cg_first,proj_axis,box_loc,the_center,width) 
    if 0: 
        print('box position: ',box)
        sigma_calc = calculate.Calculate(field_main,the_region,cg_first,cg_second) 
    if 1:
        mom_map = momentmap.moment_map(cg_first)

    if 0:
        # the following could be split into two definitions/files PDF and FIT
        bins, cbins, pdf = make_pdfNfit.Make_PDFnFIT(field_main,projection,proj_axis,cg_first,cg_second,cg_third) 
        sigma_sq = make_sigmasquared.Make_SigmaSquared(field_main,bins,cbins,pdf)
        if projection == False: 
            sigma_sq = None
    if 0:
        four_amps, powerspec2D, powerspec3D, infxyz, ks = powerspec.PowerSpec(cg_first,cg_second,dimension)

    # SIGNAL TO NOISE RATIOS
    if 0:
        collect = defaultdict(list)
        for target in target_snr:
            for many in range(100):  
                if target != inf:
                    noisy_powerspec2D, noisy_powerspec3D, noise_sign, noise_signs = signaltonoise.signaltonoise(cg_first, dimension, target, four_amps, powerspec2D, powerspec3D) 
                    #collect['noise'].append(noise_sign)
                    #collect['signal_noise'].append(noise_signs)
                if 0:
                    if target == inf:
                        noise = False
                    sigma_rho = parsevals.parsevals_th(dimension,powerspec2D,powerspec3D,infxyz,ks,sigma_sq,noise,noisy_powerspec2D,noisy_powerspec3D) 
                    collect[target].append(sigma_rho)
            
            #pdb.set_trace()
            if 0:  
                # plot distributions per SNR 
                length = len(collect[target])
                the_bins = int(np.sqrt(length))

                #the_PDFN, xbinsN = np.histogram(collect['noise'], bins=the_bins, weights = None, density=True)
                #bin_centersN = 0.5*(xbinsN[1:]+xbinsN[:-1])
                #the_PDFSNR, xbinsSRN = np.histogram(collect['signal_noise'], bins=the_bins, weights = None, density=True)
                #bin_centersSNR = 0.5*(xbinsSRN[1:]+xbinsSRN[:-1])
                the_PDFSR, xbinsSR = np.histogram(collect[target], bins=16, weights = None, density=True)
                bin_centersSR = 0.5*(xbinsSR[1:]+xbinsSR[:-1])
                # PDF(N)
                #plt.plot(bin_centersN,the_PDFN, c='orange')  
                # PDF(SNR)
                #plt.plot(bin_centersSNR,the_PDFSNR, c='blue')
                # PDF(sigma_rhos)
                plt.plot(bin_centersSR,the_PDFSR, c='blue')
                # labels
                #plt.xlim(-1,1)
                #plt.xscale('log')
                #plt.xlabel('N: orange, SNR: blue')
                plt.xlabel('sigma rho')
                plt.ylabel('PDF')
                plt.savefig('PDF_sigmasrho')

        # CALLING FOR A BOXPLOT
        if 0:
            #fig, ax=plt.subplots(1,1)
            labels, data = collect.keys(), collect.values()
            plt.boxplot(data,showmeans=True)
            '''
            if boxes[num] == 'l_box':
                plt.plot([1,6],[2.639,2.639],'--',c=[0.5]*4)
            if boxes[num] == 'c_box':
                plt.plot([1,6],[0.330,0.330],'--',c=[0.5]*4)
            if boxes[num] == 'r_box':
                plt.plot([1,6],[2.637,2.637],'--',c=[0.5]*4)
            '''
            plt.xlabel('SNR')
            plt.ylabel(r'$\sigma_{\rho/\rho_o} $')
            plt.xticks(range(1, len(labels) + 1), labels)
            #plt.ylim(0.32,0.475)
            plt.savefig('boxplot_100n_snr_%s'%boxes[num])
            plt.clf()
            # sample from kappa story: SEE MORE IN .PY FILE
            # boxes = ax1.boxplot(data, meanprops=meanpointprops, showmeans=True, meanline=False, showfliers=False)
        
        # THE PLOT
        if 0:
            print('plotting')
            #target_snr_sort = np.sort(target_snr)  #but then I'll need the sigmas array to be sorted accordingly
            #sigmas_sort = [x for _,x in sorted(zip(target_snr,sigmas))]
            #target_sort = np.sort(target_snr)
            #fig, ax = plt.subplots(1,1)
            ax.scatter(target_snr, sigmas,c=colors[num],alpha=0.7)
            ax.set(xscale='log', xlim=[min(target_snr),max(target_snr)])
            #pdb.set_trace()
            limit[0] = min([limit[0],min(target_snr)])
            limit[1] = max([limit[1],max(target_snr)])
            limit[2] = min([limit[2],min(sigmas)])
            limit[3] = max([limit[3],max(sigmas)])
            #ax.scatter(target_sort, sigmas_sort,c=colors[num],alpha=0.7)
            #print('snr_db',target_snr_db)
            #print('snr',target_snr)
            #print('sigmas',sigmas)
            
# OVERLAY THE THREE BOXES
if 0:
    ax.set_xlabel(r'$SNR$')
    ax.set_xscale('log')
    ax.set(xlim=[limit[0],limit[1]],ylim=[limit[2],limit[3]])
    ax.set_xlim(0.5,110)
    ax.set_ylim(-0.25,6.5)
    ax.set_ylabel(r'$\sigma_{\rho Brunt}$')
    outname = 'sigma_snr_fullcubes_trimbox_log_revisedfixnot'
    fig.savefig(outname)
    #fig.savefig('sigma_snr_fullcubes_log')        



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
    obtain_b.Obtain_b()  #sigma_rho,sigma_v,polytropic_index)

