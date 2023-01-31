
from startlaser import *
'''
take the fourier transform of the 2D field to get:
P2D --> P3D --> 3D variance

Note: use rfft and rfftfreq to calculate positive side only
'''

def parsevals_th(dim, p2d, p3d, dk, sigmaN_sq=None):

    if dim == 1:
        P2D = (p2d * dk).sum()
        P3D = (p3d * dk).sum()

        the_R = P2D/P3D
        sigmaRho_sq = sigmaN_sq * (P3D/P2D)
        sigmaRho = np.sqrt(sigmaRho_sq)

    if dim == 2:
        P3D = (p3d * dk).sum()
        sigmaRho = np.sqrt(P3D)
        
    #pdb.set_trace()
    if 1:
        alphaFile = open("sflaser_records.txt",'a')
        #alphaFile.write("\nThe R value %f"%(the_R))
        alphaFile.write("\nSigma 3D from P3D %f"%(sigmaRho))
        #alphaFile.write("\nSigma 3D_2D %f"%(sigmaRho))
        alphaFile.close()

    return 
