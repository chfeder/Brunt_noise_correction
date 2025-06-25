
'''
adding signal to the image (as opposed to the DFT)
'''
from startlaser import *
import fake_powerlaw    
reload(fake_powerlaw)


def add_noise(cg_main, the_snr, noisetype=None):

    if noisetype == "one": 
        signal = cg_main.v
        mean_noise = 0 
        # prep noise for the density
        scale = np.std(signal)/the_snr   #was scale = signal/snr  
        # this adds gaussian noise with (mean, standard devation, and shape)
        noise = np.random.normal(mean_noise, scale=scale, size=signal.shape)  # will be different every time, so we can run it multiple times to get a PDF(N) 
        sn_signal = signal + noise  
        if the_snr == 1000000:
            sn_signal = signal
        if 0:
            plt.clf()
            plt.hist(noise.flatten(),bins=100)
            plt.savefig('noise_kone_gauss.png')

    if noisetype == "half":  
        the_signal = cg_main.v
        theN = the_signal.shape[0]
        theK = 2.5  #1 for k=0, 1.5 for k=0.5, 2 for k=1, 2.5 for k=1.5, 3 for k=2, 4 for k=3 
        theMin = 0
        theMax = 51
        theAmp = 1
        noise = fake_powerlaw.fake_powerlaw(theN, theK, theMin, theMax, Amplitude=theAmp, phase=True)
        scale = np.std(the_signal)/the_snr   
        scale /= noise.std()  
        noise *= scale      
        sn_signal = the_signal + noise  
        if 0:
            plt.clf()
            plt.hist(noise.flatten(),bins=100)
            plt.savefig('noise_kone_fake.png')

    return sn_signal, noise 
