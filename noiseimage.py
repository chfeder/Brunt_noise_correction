
'''
adding signal to the image (as opposed to the DFT)
'''
from startlaser import *


def add_noise(cg_main, the_snr):
    signal = cg_main.v

    if the_snr != inf: 
        mean_noise = 0
        
        # prep noise for the density
        scale = np.std(signal)/the_snr   #was scale = signal/snr 
        # this adds gaussian noise with (mean, standard devation, and shape)
        noise = np.random.normal(mean_noise, scale=scale, size=signal.shape)  # will be different every time, so we can run it multiple times to get a PDF(N) 
        sn_signal = signal + noise  
    else:
        print('with snr 1000000')
        sn_signal = signal
        mean_noise = 0
        scale = np.std(signal)/1000000   #was scale = signal/snr 
        noise = np.random.normal(mean_noise, scale=scale, size=signal.shape)  # will be different every time, so we can run it multiple times to get a PDF(N) 
        #noise = 0

    return sn_signal, noise 


def imagenoise_old(cg_main, cg_two, the_snr, box_loc, imaging): 
    # a common shape: 102, 102
    if imaging == 'yes':
        if the_snr == inf:  #image the original image 
            plt.imshow(cg_main.v,interpolation='nearest',origin='lower')
            plt.savefig('imgRho_proj0375_%s.png'%(box_loc))
            plt.close()
        if the_snr != inf:  #image the noise 
            plt.imshow(noise,interpolation='nearest',origin='lower')
            plt.savefig('noise_proj0375_%d_%s.png'%(the_snr,box_loc))
            plt.close()
        if the_snr != inf:  #image the new image 
            plt.imshow(sn_signal,interpolation='nearest',origin='lower')
            plt.savefig('imgRhoNoise_proj0375_%d_%s.png'%(the_snr,box_loc))
            plt.close()
                                
    return 
