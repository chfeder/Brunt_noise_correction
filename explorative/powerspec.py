
'''
insights from the cloud tutorial: https://bertvandenbroucke.netlify.app/2019/05/24/computing-a-power-spectrum-in-python/
take the fourier transform of the 2D field to get:
P2D --> P3D --> 3D variance
Note: use rfft and rfftfreq to calculate positive side only
'''
from startlaser import *


def powerspec_lowpass(cg_main, justnoise=False, the_k=0, lowpass_filter='no'):

    npix = cg_main.shape[0] #102 

    fourier_image = np.fft.fftn(cg_main)
    fourier_amps = np.abs(fourier_image)**2
    fourier_amplitudes = fourier_amps.flatten()
 
    kfreq = np.fft.fftfreq(npix) * npix
    kfreq2D = np.meshgrid(kfreq, kfreq)  # to convert into two dim array
    knrm2Ds = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)  # equiv to np.hypot(ix,iy) 

    knrm2D = knrm2Ds.flatten()
    kbins = np.arange(0,knrm2D.max()+1,1)  
    kvals = 0.5 * (kbins[1:] + kbins[:-1])

    # the actual low pass filter
    lowpass = fourier_image.copy()
    lowpass[knrm2Ds > the_k] = 0.0
    lowpass_mag = np.abs(lowpass)**2
    fourier_lowpassamps = lowpass_mag.flatten()

    # passign the baton!!
    if lowpass_filter == 'yes':
        fourier_amplitudes = fourier_lowpassamps
 
    Abins2D, _, _ = stats.binned_statistic(knrm2D, fourier_amplitudes, \
                                        statistic = "sum", bins = kbins)  
    powerspec2D = Abins2D 
    powerspec3D = 2 * kvals * Abins2D     
    dxyz = 1/fourier_image.size 

    return kvals, powerspec2D, powerspec3D, dxyz


def powerspec(cg_main, dim, snr=0, iters=1, justnoise=False):

    # pixel resolution of the covergrid data array
    npix = cg_main.shape[0] 
    if justnoise == False:
        cg_main = cg_main/cg_main.mean() - 1   
    
    ''' DEBUG
    plt.clf()
    plt.imshow(cg_main)
    plt.savefig('density_s1')
    '''

    # the fourier transform of the covergrid data array
    if dim == 1: 
        fourier_image = np.fft.fftn(cg_main)
    if dim == 2:
        fourier_image = np.fft.fftn(cg_main) 
        #fourier_image = np.fft.fftn(cg_main,(307,307,307))  #attempt: pads the y direction with zeros

    # assuming the average amplitude is zero, which is not required, zero-th term does not show up anyway...
    # to get the size of the amplitudes only
    fourier_amps = np.abs(fourier_image)**2

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
    ''' DEBUG
    plt.clf()
    plt.imshow(fourier_amps.real)
    plt.savefig('density_fourier_s1')
    ''' 
    fourier_amplitudes = fourier_amps.flatten()
    #corresponding k values are the midpoints of the bins
    kvals = 0.5 * (kbins[1:] + kbins[:-1])
    
    # for 2D 
    if dim == 1:
        # the summed fourier amplitude squared (power) in each bin; 
        Abins2D, _, _ = stats.binned_statistic(knrm2D, fourier_amplitudes, \
                                            statistic = "sum", bins = kbins)  #this used to be "mean", now "sum"
        # total variance, aka total power
        #Abins2D *= np.pi * (kbins[1:]**2 - kbins[:-1]**2)  
        powerspec2D = Abins2D 
        powerspec3D = 2 * kvals * Abins2D 
    # for 3D
    if dim == 2:
        Abins3D, _, _ = stats.binned_statistic(knrm3D, fourier_amplitudes, \
                                            statistic = "sum", bins = kbins)  #this used to be "mean", now "sum"
        # total variance, aka total power
        #Abins3D = Abins3D * (4./3.) * np.pi * (kbins[1:]**3 - kbins[:-1]**3)   
        powerspec3D = Abins3D  
        powerspec2D = 0
    
    # define dk
    dxyz = 1/fourier_image.size 


    # checkpoints from what we've gathered; early valuable tests
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




# the following were only tests for understanding; they do not get used.
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
    # THIS WAS JUST A TEST

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