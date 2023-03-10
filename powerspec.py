
from startlaser import *
'''
insights from the cloud tutorial
'''

def PowerSpec(cg_main,cg_two,dim):
    if 0:
        the_region = mpimg.imread("NIF_hdf5_plt_cnt_0375_Slice_z_density_region.png")  
        #shape for image: (1186, 1085, 4)
        #shape for cg_main: (307, 89, 307), shape for 'square' analysis region: (179, 89, 179), 
        #actual square: (89,89,89) 
        the_region2D = mpimg.imread("imgRho_proj0375_1.png") 
        #shape for image: (480, 640, 4)
        #shape for cg_main: (307, 307), shape for square region: ()


    # PIXEL RESOLUTION OF THE COVERGRID DATA ARRAY
    npix = cg_main.shape[0] 

    # THE FOURIER TRANSFORM OF THE COVERGRID DATA ARRAY
    if dim == 1: 
        fourier_image = np.fft.fftn(cg_main)
    if dim == 2:
        fourier_image = np.fft.fftn(cg_main) 
        #fourier_image = np.fft.fftn(cg_main,(307,307,307))  #attempt: pads the y direction with zeros

    # ASSUMING THE AVERAGE AMPLITUDE IS ZERO; not required; zero-th term does not show up anyway...
    # TO GET THE SIZE OF THE AMPLITUDES ONLY
    fourier_amps = np.abs(fourier_image)**2

    # RETURNS ONE DIM ARRAY CONTAINING THE WAVE VECTORS OF THE FT CALL, IN THE CORRECT ORDER
    # * NPIX CONVERTS IT TO PIXEL FREQUENCY  
    kfreq = np.fft.fftfreq(npix) * npix
    if dim == 1:
        # TO CONVERT INTO TWO DIM ARRAY
        kfreq2D = np.meshgrid(kfreq, kfreq)
        # TO GET THE NORM
        knrm2Ds = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)
        knrm2D = knrm2Ds.flatten()
        kbins = np.arange(0,knrm2D.max()+1,1)  #2nd bins edit
    if dim == 2:
        # TO CONVERT INTO THREE DIM ARRAY
        kfreq3D = np.meshgrid(kfreq, kfreq, kfreq)
        # TO GET THE NORM
        knrm3D = np.sqrt(kfreq3D[0]**2 + kfreq3D[1]**2 + kfreq3D[2]**2)
        knrm3D = knrm3D.flatten()
        kbins = np.arange(0,knrm3D.max()+1,1)  #2nd bins edit
    
    fourier_amplitudes = fourier_amps.flatten()
    # 'CONSIDERING THE MAX POSSIBLE K VALUES POSES A SAMPLING PROBLEM': but we're still will choose that above
    #kbins = np.arange(0, npix//2+1, 1.)  #1st bins edit 
    #kbins = np.arange(0.5, npix//2+1, 1.) #bins as in the article

    # CORRESPONDING K VALUES ARE THE MIDPOINTS OF THE BINS
    kvals = 0.5 * (kbins[1:] + kbins[:-1])

    if dim == 1:
        # AVERAGE FOURIER AMPLITUDE SQUARED (POWER) IN EACH BIN
        Abins2D, _, _ = stats.binned_statistic(knrm2D, fourier_amplitudes, \
                                            statistic = "sum", bins = kbins)  #this used to be "mean", now "sum"
        # TOTAL VARIANCE, AKA TOTAL POWER 
        #Abins2D *= np.pi * (kbins[1:]**2 - kbins[:-1]**2)  #turned off
        powerspec2D = Abins2D  #the 2pi doesn't make a difference in the Brunt method 
        powerspec3D = 2 * kvals * powerspec2D 
    if dim == 2:
        # AVERAGE FOURIER AMPLITUDE SQUARED (POWER) IN EACH BIN
        Abins3D, _, _ = stats.binned_statistic(knrm3D, fourier_amplitudes, \
                                            statistic = "sum", bins = kbins)  #this used to be "mean", now "sum"
        # TOTAL VARIANCE, AKA TOTAL POWER 
        #Abins3D = Abins3D * (4./3.) * np.pi * (kbins[1:]**3 - kbins[:-1]**3)   #turned off
        powerspec3D = Abins3D  
        powerspec2D = 0
    
    # DEFINE DK 
    dxyz = 1/fourier_image.size 
    # PLANCHEREL
    if 0:
        cg_sq = (cg_main **2).sum() 
        fft_sq = (fourier_amplitudes * dxyz).sum()  
        powerspec = (powerspec3D * dxyz).sum()
        # FFT_SQ AND POWERSPEC ARE 22442.756611192683 and 22442.756611192686 respectively

        ratio = fft_sq/cg_sq  #yep, it are 1
        ratio_2 = powerspec/cg_sq  #yep, it is 1 


        # PLAY
        cg_mean = cg_main.mean().v
        cg_mean_cv = ( (cg_main * cg_two).sum()/cg_two.sum() ).v

        powerspec_sq = np.sqrt(powerspec3D)
        sigma_rho_i = (powerspec_sq * dxyz).sum()  #0.421512 
        sigma_rho_ii = np.sqrt(((powerspec_sq * dxyz).sum()**2) - cg_mean**2)  #0.3851401 
        sigma_rho_iii = (powerspec_sq * dxyz).sum() - cg_mean  #mathematically doesn't make sense, 0.250225   

        sigma_rho_iv = np.sqrt(powerspec - cg_mean**2)  #149, cg_mean**2 here is insignificant 
        print(sigma_rho_i)
        
        '''
        # exactly what IS in calculate.py...........        
        cg_mean_cv = ( (cg_main * cg_two).sum()/cg_two.sum() )
        cg_main_s = cg_main/cg_mean_cv
        cg_varsq = (((cg_main_s - 1)**2) * cg_two).sum()/cg_two.sum()  #not the same as not doing - 1
        cg_var = np.sqrt(cg_varsq)
        '''


    # PLOTTINGS SETUP
    if 0: 
        fig,ax = plt.subplots(1,1)
        fig.subplots_adjust(wspace=0, hspace=0)
                
        cg_unitless = cg_main.v
        #ax[0].imshow(cg_unitless)
        #ax.loglog(kvals,powerspec2D)
        ax.loglog(kvals,powerspec3D,c='green')
        plt.xlabel('$k$')
        plt.ylabel('$P(k)$')
        plt.tight_layout()
        plt.savefig('sort_plots/boxregion_powerspec3D.png')
 
    return powerspec2D, powerspec3D, dxyz, kvals, fourier_amps, knrm2Ds  
