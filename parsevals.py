
from startlaser import *
import turbulence_driving_analysis
reload(turbulence_driving_analysis)


def parsevals_th(dim, p2d, p3d, dk, sigmaN_sq=None):

    if dim == 1:
        if 1:
            P2D = ((p2d[1:]*dk).sum()) *dk
            P3D = ((p3d[1:]*dk).sum()) *dk
            sigmaRho_sq = sigmaN_sq * (P3D/P2D)
            sigmaRho = np.sqrt(sigmaRho_sq)

            the_R = P2D/P3D  
            the_R_sqrt = np.sqrt(P2D/P3D)

            print('sigma 2 from Brunts', sigmaRho_sq)
            print('sigma from Brunts', sigmaRho)  #orig region: 0.15031, box region: 0.108819
                                         #with cg/cgmean -1, orig: 0.33404, box: 0.495
            print('the_R', the_R) 
                                   
            print('the_R^sqrt', the_R_sqrt) 


        # testing the turbulence_driving_analysis script, 
        if 0:
            CF_sqrt_R = turbulence_driving_analysis.get_Brunt_sqrt_R(psm,km)  
            sigmaRho = np.sqrt(sigmaN_sq)/CF_sqrt_R   

        if 0:
            print('the_sqrt_R',CF_sqrt_R)  #original: 2.08746
            print('sigmaRho',sigmaRho)  #original: 0.0651347

    if dim == 2: 
        P3D = ((p3d[1:] * dk).sum()) *dk   #the same if [1:] is taken out
        sigmaRho = np.sqrt(P3D)
        print('P3D', P3D)
        print('sigma from P3D', sigmaRho)

        
    if 0:  # write out what you wish
        alphaFile = open("sflaser_records.txt",'a')
        alphaFile.write("\nSigma 3D from P3D %f"%(sigmaRho))
        alphaFile.close()

    return 


# a totally kidding attempt:
'''
here it looks as if the_R is P2D/P3D and the_R^(1/2) is np.sqrt(P2D/P3D)
and Christoph "typo" comment does not make sense

# new idea
if 0:  #BUT THIS DOES NOT COMMUTE, NVM
    P2D = (np.sqrt(p2d)*dk).sum()
    P3D = (np.sqrt(p3d)*dk).sum()
    sigmaRho = np.sqrt(sigmaN_sq) * (P3D/P2D)  # box: 0.211, original: 0.435
    
    # calculated this way the_R is P2D/P3D, each Ps already square rooted, the 1/2 denoting how it was computed.
    # and Christoph's "typo" comment makes sense 
    the_R = P2D/P3D  # box:0.49076728902 , original: 0.31230609
'''