
'''
adding signal to the image (as opposed to the DFT)
'''
from startlaser import *


def add_noise(cg_main, the_snr):
    signal = cg_main.v

    mean_noise = 0
    
    # prep noise for the density
    scale = np.std(signal)/the_snr   #was scale = signal/snr 
    
    # this adds gaussian noise with (mean, standard devation, and shape)
    noise = np.random.normal(mean_noise, scale=scale, size=signal.shape)  # will be different every time, so we can run it multiple times to get a PDF(N) 

    sn_signal = signal + noise  

    if the_snr == 1000000:
        sn_signal = signal

    return sn_signal, noise 
