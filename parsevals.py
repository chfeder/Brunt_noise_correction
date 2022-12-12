
'''
take the fourier transform of the 2D field to get:
P2D --> P3D --> 3D variance

Note: use rfft and rfftfreq to calculate positive side only
'''
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import fft
#from scipy import fftfreq  #combining scipy and np for now...hm


def parsevals_th(field2D,sigmaN_sq):

    n = field2D.size 
    print('field size',n)
    p2d = fft(field2D)
    print('p2d',p2d)
    p2dx = np.fft.fftfreq(n)   #could add ,d=sample spacing
    print('p2dx',p2dx)
    wavenumber = 1/p2dx 
    print('wavenumber', wavenumber)

    p3d = 2 * wavenumber * p2d
    print('p3d',p3d)
    sigmaRho_sq = (sigmaN_sq * p3d)/p2d

    if 1:
        alphaFile = open("sflaser_records.txt",'a')
        alphaFile.write("\nSigma^2 3D %f"%(sigmasquared))
        alphaFile.close()

    return 