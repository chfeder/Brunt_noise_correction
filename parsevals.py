
'''
take the fourier transform of the 2D field to get:
P2D --> P3D --> 3D variance

Note: use rfft and rfftfreq to calculate positive side only
'''
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import fft
import pdb
#from scipy import fftfreq  #combining scipy and np for now...hm


def parsevals_th(dim,sigmaN_sq, p2d, p3d, dk):

    if dim == 1:
        P2D = (p2d * dk).sum()
        P3D = (p3d * dk).sum()
        sigma_rho_P3D = np.sqrt(P3D) 
        print('sigma_rho_P3D',sigma_rho_P3D)
    if dim == 2:
        p2d = 64130704.256487735
        P2D = p2d
        P3D = (p3d * dk).sum()
        sigma_rho_P3D = np.sqrt(P3D) 
        print('sigma_rho_P3D',sigma_rho_P3D)

    sigmaRho_sq = sigmaN_sq * (P3D/P2D)
    sigmaRho = np.sqrt(sigmaRho_sq)
    #pdb.set_trace()

    if 0:
        alphaFile = open("sflaser_records.txt",'a')
        alphaFile.write("\nP2D: %f, P3D_2D: %f, P3D: %f"%(P2D,P3D,P3D))
        alphaFile.close()
    if 0:
        alphaFile = open("sflaser_records.txt",'a')
        alphaFile.write("\nSigma 3D %f"%(sigmaRho))
        alphaFile.close()

    return 