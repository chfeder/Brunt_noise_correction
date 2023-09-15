
'''
adding signal to the image (as opposed to the FT)
'''

from startlaser import *


def noiseimage(cg_main, cg_two, the_snr, box_loc, imaging): 

    # shape 102, 102
    signal = cg_main.v
    signal_vol = cg_two.v

    if the_snr != inf: 
        mean_noise = 0
        
        # prep noise to the density
        scale = np.std(signal)/the_snr   #was scale = signal/snr 
        noise = np.random.normal(mean_noise, scale=scale, size=signal.shape)  # will be different every time, so we can run it multiple times to get a PDF(N) 
        sn_signal = signal + noise  
        #stop()
        '''
        # prep noise to the volume
        scale_two = np.std(signal_vol)/the_snr
        noise_two = np.random.normal(mean_noise, scale=scale_two, size=signal_vol.shape)  
        sn_vol = signal_vol + noise
        '''
        sn_vol = signal_vol

    else:
        print('this should not be printed')
        sn_signal = signal
        sn_vol = signal_vol 
        noise = 0


    if 0:  #update on imaging 
        if imaging == 'yes':
            if the_snr == inf:  #image the original image 
                cfp.plot_map(signal,save='signal_inf_%s.png'%box_loc)
            if the_snr != inf:  #image the noise & noisy image
                cfp.plot_map(noise,save='noise_%d_%s.png'%(the_snr,box_loc))
                cfp.plot_map(sn_signal,save='noiseimage_%d_%s.png'%(the_snr,box_loc))
    if 0:
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
                                
    return sn_signal, sn_vol, noise 
