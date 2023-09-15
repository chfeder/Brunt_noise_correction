
'''
explore signal to noise ratios
'''

from startlaser import *


def signaltonoise(cg_main, dim, the_snr, signal, power2D=0, power3D=0): 

    if dim == 1:  #to be edited for 3D
        powerspec = power2D
    if dim == 2:
        powerspec = power3D 
    if 1: 
        #np.random.seed(4545118)  #to keep it fixed; must be between 0 and 2**32-1
        mean_noise = 0
        scale = np.std(signal)/the_snr  #was scale = signal/the_snr
        noise = np.random.normal(mean_noise, scale=scale, size=signal.shape)  # will be different every time, so we can run it multiple times to get a PDF(N)
                                                                             # with sigma=S/SNR ... 
        sn_signal = signal + noise  # on tutorial, they add it to the amplitude signal...but then the shape must be of the signal; Dr. Collins agreed
                                    # ...so then we will also have a PDF(SN), msg says PDF(S)...shifted but with mean value sitting on original S


    # FROM powerspec.py  
    npix = cg_main.shape[0]
    kfreq = np.fft.fftfreq(npix) * npix
    if dim ==1:
        kfreq2D = np.meshgrid(kfreq, kfreq)  
        knrm2Ds = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)  
        knrm2D = knrm2Ds.flatten()
        kbins = np.arange(0,knrm2D.max()+1,1)   
        kvals = 0.5 * (kbins[1:] + kbins[:-1])

        snrpower2D, _, _ = stats.binned_statistic(knrm2D, sn_signal, \
                                            statistic = "sum", bins = kbins)  #this used to be "mean", now "sum"
        snrpower3D = 2 * kvals * snrpower2D 
    if dim == 2:
        kfreq3D = np.meshgrid(kfreq, kfreq, kfreq)
        knrm3D = np.sqrt(kfreq3D[0]**2 + kfreq3D[1]**2 + kfreq3D[2]**2)
        knrm3D = knrm3D.flatten()
        kbins = np.arange(0,knrm3D.max()+1,1) 

        snrpower3D, _, _ = stats.binned_statistic(knrm3D, sn_signal, \
                                            statistic = "sum", bins = kbins)  
        snrpower2D = 0



    #PLAY
    if 0:  #pure signal with signal plus noise
        fig,ax = plt.subplots(1,1)
        fig.subplots_adjust(wspace=0, hspace=0)
                
        cg_unitless = cg_main.v 
        kvals = 0.5 * (kbins[1:] + kbins[:-1])

        ax.loglog(kvals,power3D,c='green')
        ax.loglog(kvals,snrpower3D,c='red')
        plt.xlabel('$k$')
        plt.ylabel('$P(k)$')
        plt.tight_layout()
        plt.savefig('p3D_snrp3D_%d.png'%the_snr) 
    if 0:  #image the noise!
        fourier_image = np.fft.fftn(cg_main)
        fourier_ampss = np.abs(fourier_image)**2
        fourier_ampss_flat = fourier_ampss.flatten()
        #pdb.set_trace()
        fourier_image_noise = fourier_ampss_flat + noise  
        fourier_image_noise_inverse = np.fft.ifft(fourier_image_noise)
        plt.imshow(fourier_image_noise_inverse.sum(0))  #along which axis? 
    if 0:  #measure snr
        def signaltonoise(a, axis=0, ddof=0):
            a = np.asanyarray(a)
            m = a.mean(axis)
            sd = a.std(axis=axis, ddof=ddof)
            return np.where(sd == 0, 0, m/sd)

        arr = snrpower3D  #sn_signal 
        snr = signaltonoise(arr)
        print('what the snr mean/std def gives with %d'%the_snr,snr)
        #pdb.set_trace()


    return snrpower2D, snrpower3D, noise, sn_signal
