
'''
explore signal to noise ratios
'''

from startlaser import *


def signaltonoise(cg_main, dim, snr_db, signal, power2D=0, power3D=0): 

    if dim == 1:  #to be edited for 3D
        powerspec = power2D
    if dim == 2:
        powerspec = power3D 

    if 1:  #snr in decibels
        sig_avg_watts = np.mean(powerspec)  #we already have our signal power, still need mean?
        sig_avg_db = 10 * np.log10(sig_avg_watts)

        noise_db = sig_avg_db - snr_db 
        noise_watts = 10 ** (noise_db/ 10)  #triple check
        
        mean_noise = 0
        np.random.seed(4545118)  #to keep it fixed; must be between 0 and 2**32-1
        noise = np.random.normal(mean_noise, np.sqrt(noise_watts), len(signal))
        sn_signal = signal + noise  #on tutorial, they add it to the amplitude signal...but then the shape must be of the signal; Dr. Collins agreed


    # from powerspec.py  
    npix = cg_main.shape[0]
    kfreq = np.fft.fftfreq(npix) * npix
    if dim ==1:
        kfreq2D = np.meshgrid(kfreq, kfreq)  
        knrm2Ds = np.sqrt(kfreq2D[0]**2 + kfreq2D[1]**2)  
        knrm2D = knrm2Ds.flatten()
        kbins = np.arange(0,knrm2D.max()+1,1)   

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



    if 0:  #pure signal with signal plus noise
        fig,ax = plt.subplots(1,1)
        fig.subplots_adjust(wspace=0, hspace=0)
                
        cg_unitless = cg_main.v 
        kvals = 0.5 * (kbins[1:] + kbins[:-1])

        ax.loglog(kvals,powerspec,c='green')
        ax.loglog(kvals,snrpower3D,c='red')
        plt.xlabel('$k$')
        plt.ylabel('$P(k)$')
        plt.tight_layout()
        plt.savefig('p3D_snrp3D_%d.png'%snr)

    return snrpower2D, snrpower3D
