
'''
adding signal to the image (as opposed to the FT)
'''

from startlaser import *


def noiseimage(cg_main, the_snr, box_loc, imaging): 

    # shape 102, 102
    #pdb.set_trace()
    signal = cg_main.v

    if the_snr != inf: 
        mean_noise = 0
        scale = signal/the_snr  #this is sigma = S/SNR
        noise = np.random.normal(mean_noise, scale=scale, size=signal.shape)  # will be different every time, so we can run it multiple times to get a PDF(N) 
        sn_signal = signal + noise  
        #pdb.set_trace()
    else:
        sn_signal = signal

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
                                
    return sn_signal 
