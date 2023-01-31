
from startlaser import *
'''
insights from the cloud tutorial
'''

def PowerSpec(cg_main, dim):
    if 0:
        the_region = mpimg.imread("NIF_hdf5_plt_cnt_0375_Slice_z_density_region.png")  
        #shape for image: (1186, 1085, 4)
        #shape for cg_main: (307, 89, 307), shape for 'square' analysis region: (179, 89, 179), 
        #actual square: (89,89,89) 
        the_region2D = mpimg.imread("imgRho_proj0375_1.png") 
        #shape for image: (480, 640, 4)
        #shape for cg_main: (307, 307), shape for square region: ()

    npix = cg_main.shape[0] 

    if dim == 1: 
        fourier_image = np.fft.fftn(cg_main)
    if dim == 2:
        fourier_image = np.fft.fftn(cg_main) 
        #fourier_image = np.fft.fftn(cg_main,(307,307,307))  #pads the y direction with zeros
    fourier_amplitudes = np.abs(fourier_image)**2

    kfreq = np.fft.fftfreq(npix) * npix
    if dim == 1:
        kfreq2D = np.meshgrid(kfreq, kfreq)
        knrm2D = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)
        knrm2D = knrm2D.flatten()
        kbins = np.arange(0,knrm2D.max()+1,1)  #edited
    if dim == 2:
        kfreq3D = np.meshgrid(kfreq, kfreq, kfreq)
        knrm3D = np.sqrt(kfreq3D[0]**2 + kfreq3D[1]**2 + kfreq3D[2]**2)
        knrm3D = knrm3D.flatten()
        kbins = np.arange(0,knrm3D.max()+1,1)  #edited

    
    fourier_amplitudes = fourier_amplitudes.flatten()

    #kbins = np.arange(0, npix//2+1, 1.)  #have to line up 
    #kbins = np.arange(0.5, npix//2+1, 1.)
    kvals = 0.5 * (kbins[1:] + kbins[:-1])

    if dim == 1:
        Abins2D, _, _ = stats.binned_statistic(knrm2D, fourier_amplitudes, \
                                            statistic = "sum", bins = kbins)  #this used to be "mean"
        #Abins2D *= np.pi * (kbins[1:]**2 - kbins[:-1]**2)  #turned off
        powerspec2D = Abins2D
        powerspec3D = 2 * kvals * powerspec2D 
        dxyz = 1/kvals
    if dim == 2:
        Abins3D, _, _ = stats.binned_statistic(knrm3D, fourier_amplitudes, \
                                            statistic = "sum", bins = kbins)
        #Abins3D *= (4./ 3.) * np.pi * (kbins[1:]**3 - kbins[:-1]**3)
        powerspec3D = Abins3D  #this works.
        powerspec2D = 0
        dxyz = 1/kvals


    # PLANCHEREL
    if 0:
        # see nazare laser/powerspec_test.py
        dk = 1/fourier_image.size
        fft_sq = ((np.abs(fourier_image) **2) * dk).sum() 
        cg_sq = (cg_main **2).sum()  
        ratio = fft_sq/cg_sq
        print('ratio',ratio)

    # PLOTTINGS SETUP
    if 0: 
        fig,ax = plt.subplots(1,2)
        fig.subplots_adjust(wspace=0, hspace=0)
        
        cg_unitless = cg_main.v
        ax[0].imshow(cg_unitless)
        ax[1].loglog(kvals,powerspec2D)
        ax[1].loglog(kvals,powerspec3D,c='green')
        plt.xlabel('$k$')
        plt.ylabel('$P(k)$')
        plt.tight_layout()
        plt.savefig('sort_plots/region_powerspec2D_3D.png')
 
    return powerspec2D, powerspec3D, dxyz, fourier_image 
