
'''
insights from the cloud tutorial: https://bertvandenbroucke.netlify.app/2019/05/24/computing-a-power-spectrum-in-python/
take the fourier transform of the 2D field to get:
P2D --> P3D --> 3D variance
Note: use rfft and rfftfreq to calculate positive side only
'''
from startlaser import *


def lowpassfilter(cg_main):
    # EXAMPLE METHOD ONLY 

    F = np.fft.fftn(cg_main)
    F_shift = np.fft.fftshift(F)  #so that k=0 is a the center
    F_mag = np.abs(F)  #amplitude
    F_phase = np.angle(F_shift)  #phase

    # filtering out high frequencies
    ix, iy = np.mgrid[0:F_shift.shape[0], 0:F_shift.shape[1]]
    ix -= F_shift.shape[0]//2
    iy -= F_shift.shape[1]//2

    F_filtered = F_shift.copy()
    F_filtered[np.hypot(ix, iy) > 25] = 0.0
    # ...imshow(np.log(np.abs(F_filtered)))

    # transform this back and see the result
    f_filtered = np.fft.ifft2(np.fft.ifftshift(F_filtered))
    # ...imshow(f_filtered.real, cmap=:w"gray")

    return


def powerspec_lowpass_test(cg_main, justnoise=False, the_k=0, lowpass_filter='no'):

    if justnoise == False:
        cg_main = cg_main/cg_main.mean() - 1   

    npix = cg_main.shape[0] 

    fourier_image = np.fft.fftn(cg_main)
    F_shift = np.fft.fftshift(fourier_image)  #so that k=0 is a the center
    #fourier_amps = np.abs(fourier_image)**2
    fourier_amps = np.abs(F_shift)**2
    fourier_amplitudes = fourier_amps.flatten()

    ix, iy = np.mgrid[0:F_shift.shape[0], 0:F_shift.shape[1]]
    ix -= F_shift.shape[0]//2
    iy -= F_shift.shape[1]//2

    F_two = np.hypot(ix, iy)
 
    '''
    kfreq = np.fft.fftfreq(npix) * npix
    kfreq2D = np.meshgrid(kfreq, kfreq)  # to convert into two dim array
    knrm2Ds = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)  # equiv to np.hypot(ix,iy) 
    '''

    #knrm2D = knrm2Ds.flatten()
    knrm2D = F_two.flatten()
    kbins = np.arange(0,knrm2D.max()+1,1)  
    kvals = 0.5 * (kbins[1:] + kbins[:-1])


    # THE ACTUAL LOW PASS FILTER
    lowpass = F_shift.copy()
    lowpass[F_two > the_k] = 0.0
    lowpass_mag = np.abs(lowpass)**2
    fourier_lowpassamps = lowpass_mag.flatten()

    # PASSING THE BATON!!
    if lowpass_filter == 'yes':
        fourier_amplitudes = fourier_lowpassamps
 
    print('went through it')
    Abins2D, _, _ = stats.binned_statistic(knrm2D, fourier_amplitudes, \
                                        statistic = "sum", bins = kbins)  
    powerspec2D = Abins2D 
    powerspec3D = 2 * kvals * Abins2D     
    dxyz = 1/fourier_image.size 

    return kvals, powerspec2D, powerspec3D, dxyz

def powerspec_lowpass(cg_main, justnoise=False, the_k=0, lowpass_filter='no'):

    #if justnoise == False:
    #    cg_main = cg_main/cg_main.mean() - 1   

    npix = cg_main.shape[0] 

    fourier_image = np.fft.fftn(cg_main)
    fourier_amps = np.abs(fourier_image)**2
    fourier_amplitudes = fourier_amps.flatten()
 
    kfreq = np.fft.fftfreq(npix) * npix
    kfreq2D = np.meshgrid(kfreq, kfreq)  # to convert into two dim array
    knrm2Ds = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)  # equiv to np.hypot(ix,iy) 

    knrm2D = knrm2Ds.flatten()
    kbins = np.arange(0,knrm2D.max()+1,1)  
    kvals = 0.5 * (kbins[1:] + kbins[:-1])


    # THE ACTUAL LOW PASS FILTER
    lowpass = fourier_image.copy()
    lowpass[knrm2Ds > the_k] = 0.0
    lowpass_mag = np.abs(lowpass)**2
    fourier_lowpassamps = lowpass_mag.flatten()

    # PASSING THE BATON!!
    if lowpass_filter == 'yes':
        fourier_amplitudes = fourier_lowpassamps
 
    Abins2D, _, _ = stats.binned_statistic(knrm2D, fourier_amplitudes, \
                                        statistic = "sum", bins = kbins)  
    powerspec2D = Abins2D 
    powerspec3D = 2 * kvals * Abins2D     
    dxyz = 1/fourier_image.size 

    return kvals, powerspec2D, powerspec3D, dxyz



def powerspec(cg_main, dim, box_loc, snr, rep):

    # pixel resolution of the covergrid data array
    npix = cg_main.shape[0] 
    cg_main = cg_main/cg_main.mean() - 1   

    # the fourier transform of the covergrid data array
    if dim == 1: 
        fourier_image = np.fft.fftn(cg_main)
        F_shift = np.fft.fftshift(fourier_image)  #so that k=0 is a the center
    if dim == 2:
        fourier_image = np.fft.fftn(cg_main) 
        #fourier_image = np.fft.fftn(cg_main,(307,307,307))  #attempt: pads the y direction with zeros

    # assuming the average amplitude is zero, which is not required, zero-th term does not show up anyway...
    # to get the size of the amplitudes only
    fourier_amps = np.abs(fourier_image)**2
    #print('using f-shift')
    #fourier_amps = np.abs(F_shift)**2

    # returns one dim array containing the wave vectors of the ft call, in the correct order
    # * npix converts it to pixel frequency 
    kfreq = np.fft.fftfreq(npix) * npix
    # for 2D
    if dim == 1:
        kfreq2D = np.meshgrid(kfreq, kfreq)  # to convert into two dim array
        knrm2Ds = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)  # to get the norm
        knrm2D = knrm2Ds.flatten()
        kbins = np.arange(0,knrm2D.max()+1,1)  
    # for 3D
    if dim == 2:  
        kfreq3D = np.meshgrid(kfreq, kfreq, kfreq)
        knrm3D = np.sqrt(kfreq3D[0]**2 + kfreq3D[1]**2 + kfreq3D[2]**2)
        knrm3D = knrm3D.flatten()
        kbins = np.arange(0,knrm3D.max()+1,1) 

    fourier_amplitudes = fourier_amps.flatten()
    '''
    ix, iy = np.mgrid[0:F_shift.shape[0], 0:F_shift.shape[1]]
    ix -= F_shift.shape[0]//2
    iy -= F_shift.shape[1]//2

    F_filtered = F_shift.copy()
    F_filtered[np.hypot(ix, iy) > 25] = 0.0
    '''
    #corresponding k values are the midpoints of the bins
    kvals = 0.5 * (kbins[1:] + kbins[:-1])
    
    # for 3D
    if dim == 1:
        # the summed fourier amplitude squared (power) in each bin; 
        Abins2D, _, _ = stats.binned_statistic(knrm2D, fourier_amplitudes, \
                                            statistic = "sum", bins = kbins)  #this used to be "mean", now "sum"
        # total variance, aka total power
        #Abins2D *= np.pi * (kbins[1:]**2 - kbins[:-1]**2)  #turned off
        powerspec2D = Abins2D 
        powerspec3D = 2 * kvals * Abins2D 
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
    #ax[0].imshow(cg_unitless)
    if 0:
        # NEXT pass this to main to bundle iterations as a single color with all SNRs in one plot!
        if rep == 0:
            print('inside rep 0')
            plt.clf() 
        fig,ax = plt.subplots(1,1)
        fig.subplots_adjust(wspace=0, hspace=0)
        '''
        limit = np.zeros(4) 
        limit[0] = min([limit[0],min(kvals)])
        limit[1] = max([limit[1],max(kvals)])
        limit[2] = min([limit[2],min(powerspec2D)])
        limit[3] = max([limit[3],max(powerspec2D)])
        ax.set(xlim=[limit[0],limit[1]],ylim=[limit[2],limit[3]])
        '''
        ax.set_yscale('log')
        ax.loglog(kvals,powerspec2D,c='blue')
        #plt.plot(kvals,powerspec2D,c='blue') 
        plt.xlabel('$k$')
        plt.ylabel('$P(k)$')
        plt.ylim(1e0,1e11)
        plt.tight_layout()
        plt.savefig('brunts_P2D_%s_%.3f_%d_loglog.png'%(box_loc,snr,rep))

    # checkpoints from what we've gathered
    if 0:
        cg_sq = (cg_main **2).sum() 
        fft_sq = (fourier_amplitudes * dxyz).sum()  
        powerspec = (powerspec3D * dxyz).sum() 

        ratio = fft_sq/cg_sq  #yep, it is 1
        ratio_2 = powerspec/cg_sq  #yep, it is 1 

        cv = cg_two
        #cg_mean = ( (cg_main*cv).sum()/cv.sum() )
        cg_mean = np.mean(cg_main)
        #sfield_ii = cg_main/cg_mean
        #sigma_sq_ii = (((sfield_ii - 1) ** 2)*cv).sum() / cv.sum()
        sfield_ii = cg_main
        sigma_sq_ii = (((sfield_ii - cg_mean) ** 2)*cv).sum() / cv.sum()
          
        P3D = ((powerspec3D[1:] * dxyz).sum()) * dxyz 
        print('P3D', P3D)

 
    return fourier_amplitudes, powerspec2D, powerspec3D, dxyz, kvals 
