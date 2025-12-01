
from startlaser import *

def parsevals_th(dim, p2d, p3d, kvalues, p2d_ave=None, p3d_ave=None, sigmaN_sq=None, the_k=100, targeting=None):

    if dim == 1:
        if the_k == 100:
            P2D = p2d[1:].sum()
            P3D = p3d[1:].sum()
        if the_k == 0:  #k=0 calculation should be dropped and this is what this does
            P2D = np.nan
            P2D_ave = np.nan
            P3D = np.nan
            P3D_ave = np.nan
        sigmaden=0
        sigmaN = 0
        if the_k != 0 and the_k!= 100:
            k = int(the_k)
            P2D = p2d[1:k].sum()
            P2D_ave = p2d_ave[1:k].sum()
            P3D = p3d[1:k].sum()
            P3D_ave = p3d_ave[1:k].sum()

            if 1:  # given that the non-normalized values are used for all but sigmaN
                sigmaN = np.sqrt(P2D_ave)
            if 0:  # given that the normalized values are used
                sigmaN = np.sqrt(P2D)
            the_R_ave = P2D_ave/P3D_ave 
            the_R_sqrt_ave = np.sqrt(the_R_ave)
            sigmaden = sigmaN * the_R_sqrt_ave**(-1) 

        sigmaRho_sq = sigmaN_sq * (P3D/P2D)
        sigmaRhoBr = np.sqrt(sigmaRho_sq)

        the_R = P2D/P3D
        the_R_sqrt = np.sqrt(the_R)

    if dim == 2:
        P3D = p3d[1:].sum()
        sigmaRho = np.sqrt(P3D)
        the_R_sqrt = 0
 
    return sigmaRhoBr, sigmaN, the_R_sqrt, sigmaden

