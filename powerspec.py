
'''
insights from the cloud tutorial

take the fourier transform of the 2D field to get:
P2D --> P3D --> 3D variance

Note: use rfft and rfftfreq to calculate positive side only
'''
from startlaser import *


def PowerSpec(cg_main,cg_two,dim):
    if 0:  # this was a test for the cloud tutorial
        the_region = mpimg.imread("NIF_hdf5_plt_cnt_0375_Slice_z_density_region.png")  
        #shape for image: (1186, 1085, 4)
        #shape for cg_main: (307, 89, 307), shape for 'square' analysis region: (179, 89, 179), 
        #actual square: (89,89,89) 
        the_region2D = mpimg.imread("imgRho_proj0375_1.png") 
        #shape for image: (480, 640, 4)
        #shape for cg_main: (307, 307), shape for square region: ()


    # pixel resolution of the covergrid data array
    npix = cg_main.shape[0] 
    cg_main = cg_main/cg_main.mean() - 1   

    # the fourier transform of the covergrid data array
    if dim == 1: 
        fourier_image = np.fft.fftn(cg_main)
    if dim == 2:
        fourier_image = np.fft.fftn(cg_main) 
        #fourier_image = np.fft.fftn(cg_main,(307,307,307))  #attempt: pads the y direction with zeros

    # assuming the average amplitude is zero, which is not required, zero-th term does not show up anyway...
    # to get the size of the amplitudes only
    fourier_amps = np.abs(fourier_image)**2

    # returns one dim array containing thee wave vectors of the ft call, in the correct order
    # * npix converts it to pixel frequency 
    kfreq = np.fft.fftfreq(npix) * npix
    # for 2D
    if dim == 1:
        kfreq2D = np.meshgrid(kfreq, kfreq)  # to convert into two dim array
        knrm2Ds = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)  # to get the norm
        knrm2D = knrm2Ds.flatten()
        kbins = np.arange(0,knrm2D.max()+1,1)   # a fixed attempt at binning 
    # for 3D
    if dim == 2:  
        kfreq3D = np.meshgrid(kfreq, kfreq, kfreq)
        knrm3D = np.sqrt(kfreq3D[0]**2 + kfreq3D[1]**2 + kfreq3D[2]**2)
        knrm3D = knrm3D.flatten()
        kbins = np.arange(0,knrm3D.max()+1,1) 

    #pdb.set_trace()
    fourier_amplitudes = fourier_amps.flatten()
    # considering the max possible k values poses a sampling problem
    #kbins = np.arange(0, npix//2+1, 1.)  #bins edit one 
    #kbins = np.arange(0.5, npix//2+1, 1.) #bins as in the article

    #corresponding k values are the midpoints of the bins
    kvals = 0.5 * (kbins[1:] + kbins[:-1])
    
    # for 3D
    if dim == 1:
        # average (now the sum) fourier amplitude squared (power) in each bin; "mean" for snr
        Abins2D, _, _ = stats.binned_statistic(knrm2D, fourier_amplitudes, \
                                            statistic = "sum", bins = kbins)  #this used to be "mean", now "sum"
        # total variance, aka total power
        #Abins2D *= np.pi * (kbins[1:]**2 - kbins[:-1]**2)  #turned off
        powerspec2D = Abins2D 
        powerspec3D = 2 * kvals * powerspec2D 
    # for 2D
    if dim == 2:
        Abins3D, _, _ = stats.binned_statistic(knrm3D, fourier_amplitudes, \
                                            statistic = "sum", bins = kbins)  #this used to be "mean", now "sum"
        #Abins3D = Abins3D * (4./3.) * np.pi * (kbins[1:]**3 - kbins[:-1]**3)   #turned off
        powerspec3D = Abins3D  
        powerspec2D = 0
    
    # define dk
    dxyz = 1/fourier_image.size 



    # let's plot
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

    # some tests on the spot from what we've gathered
    if 0:
        cg_sq = (cg_main **2).sum() 
        fft_sq = (fourier_amplitudes * dxyz).sum()  
        powerspec = (powerspec3D * dxyz).sum() 
        print('cg_sq',cg_sq*dxyz)
        print('fft_sq',fft_sq)
        print('fft_sq shape',fft_sq)
        print('powerspec',powerspec)
        print('powerspec shape',powerspec3D.shape)
        print('fourier amps',fourier_amplitudes.shape)

        ratio = fft_sq/cg_sq  #yep, it is 1
        ratio_2 = powerspec/cg_sq  #yep, it is 1 

        cv = cg_two
        #cg_mean = ( (cg_main*cv).sum()/cv.sum() )
        cg_mean = np.mean(cg_main)
        print('cg_mean',cg_mean)
        #sfield_ii = cg_main/cg_mean
        #sigma_sq_ii = (((sfield_ii - 1) ** 2)*cv).sum() / cv.sum()
        sfield_ii = cg_main
        sigma_sq_ii = (((sfield_ii - cg_mean) ** 2)*cv).sum() / cv.sum()
        print('sigma_sq_ii',sigma_sq_ii)
        print('cg_main.var',cg_main.var())
          
        #print('p3D',powerspec3D)
        P3D = ((powerspec3D[1:] * dxyz).sum()) * dxyz 
        print('P3D', P3D)

 
    return fourier_amplitudes, powerspec2D, powerspec3D, dxyz, kvals 
