# LJV
The following numbered notes come from developing the main.py file.
All modules to be imported ar in startlaser.py 
Start with main.py  


# FROM main.py
0.
THE DATA: the simulations are in /scratch/ek9/ccf100/nif
* 9 different simulations with a few snapshots:
* NIF_hdf5_plt_cnt_0*
* frames = [000,125,250,275,300,325,350,375,425]

1.
# it was considered to do a derived field for the density to take the mean out, aka the k=0 mode
def add_field(thing):
    def overdensity(field,data):
        rho_0 = data.get_field_parameter('rho_0')
        #rho_ratio = (data['density']/rho_0).v
        #return (rho_ratio-1)
        # OR SAME THING:
        return(((data['density'] - rho_0)/rho_0))
    thing.add_field(YT_overdensity,function=overdensity)#,units='dimensionless')
add_field(the_ds)

2. 
# original target signal to noise suggested ratios: inf,100(20db),10(10db), 5(6.9897db), 3(4.7712db), 1(0db)

3. 
# attempt at distributions per SNR
    #collect['noise'].append(noise_sign)
    #collect['signal_noise'].append(noise_signs)

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

4. 
# boxplot samples: 
sample from kappa story: SEE MORE IN .PY FILE
boxes = ax1.boxplot(data, meanprops=meanpointprops, showmeans=True, meanline=False, showfliers=False)

5. 
# initial attempt at plotting singular sigma rhos: "limit" lesson
    #target_snr_sort = np.sort(target_snr)  #but then I'll need the sigmas array to be sorted accordingly
    #sigmas_sort = [x for _,x in sorted(zip(target_snr,sigmas))]
    #target_sort = np.sort(target_snr)
    #fig, ax = plt.subplots(1,1)
    ax.scatter(target_snr, sigmas,c=colors[num],alpha=0.7)
    ax.set(xscale='log', xlim=[min(target_snr),max(target_snr)])
    limit[0] = min([limit[0],min(target_snr)])
    limit[1] = max([limit[1],max(target_snr)])
    limit[2] = min([limit[2],min(sigmas)])
    limit[3] = max([limit[3],max(sigmas)])
    #ax.scatter(target_sort, sigmas_sort,c=colors[num],alpha=0.7)

    ax.set_xlabel(r'$SNR$')
    ax.set_xscale('log')
    ax.set(xlim=[limit[0],limit[1]],ylim=[limit[2],limit[3]])
    ax.set_xlim(0.5,110)
    ax.set_ylim(-0.25,6.5)
    ax.set_ylabel(r'$\sigma_{\rho Brunt}$')
    outname = 'sigma_snr_fullcubes_trimbox_log_revisedfixnot'
    fig.savefig(outname)


# CLOUD TUTORIAL IDEAS:
    if 0:  # this was a test for the cloud tutorial
        the_region = mpimg.imread("NIF_hdf5_plt_cnt_0375_Slice_z_density_region.png")  
        #shape for image: (1186, 1085, 4)
        #shape for cg_main: (307, 89, 307), shape for 'square' analysis region: (179, 89, 179), 
        #actual square: (89,89,89) 
        the_region2D = mpimg.imread("imgRho_proj0375_1.png") 
        #shape for image: (480, 640, 4)
        #shape for cg_main: (307, 307), shape for square region: ()

    # considering the max possible k values poses a sampling problem
    kbins = np.arange(0, npix//2+1, 1.)  #bins edit one 
    kbins = np.arange(0.5, npix//2+1, 1.) #bins as in the article