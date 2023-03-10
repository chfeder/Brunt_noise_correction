
from startlaser import *
import turbulence_driving_analysis
reload(turbulence_driving_analysis)
'''
take the fourier transform of the 2D field to get:
P2D --> P3D --> 3D variance

Note: use rfft and rfftfreq to calculate positive side only
'''

def parsevals_th(dim, p2d, p3d, dk,k, psm, km, sigmaN_sq=None):

    if dim == 1:
        # originally
        if 1:
            P2D = (p2d*dk).sum()
            P3D = (p3d*dk).sum()
            sigmaRho_sq = sigmaN_sq * (P3D/P2D)
            sigmaRho = np.sqrt(sigmaRho_sq)
            
            # here it looks as if the_R is P2D/P3D and the_R^(1/2) is np.sqrt(P2D/P3D)
            # and Christoph "typo" comment does not make sense
            the_R_sqrt = np.sqrt(P2D/P3D)

            print('sigmaRho', sigmaRho)  #orig region: 0.15031, box region: 0.108819
            print('the_R_sqrt', the_R_sqrt)  #orig region: 0.90455, box region: 0.95357421
        
        # new idea
        if 0:  #BUT THIS DOES NOT COMMUTE, NVM
            P2D = (np.sqrt(p2d)*dk).sum()
            P3D = (np.sqrt(p3d)*dk).sum()
            sigmaRho = np.sqrt(sigmaN_sq) * (P3D/P2D)  # box: 0.211, original: 0.435
            
            # calculated this way the_R is P2D/P3D, each Ps already square rooted, the 1/2 denoting how it was computed.
            # and Christoph's "typo" comment makes sense 
            the_R = P2D/P3D  # box:0.49076728902 , original: 0.31230609
            print('the_R',the_R)
        
        # turbulence_driving_analysis, 
        if 1:
            CF_sqrt_R = turbulence_driving_analysis.get_Brunt_sqrt_R(psm,km)  
            sigmaRho = np.sqrt(sigmaN_sq)/CF_sqrt_R   

        if 1:
            print('the_sqrt_R',CF_sqrt_R)  #original: 2.08746
            print('sigmaRho',sigmaRho)  #original: 0.0651347

    if dim == 2: 
        # SEE POWERSPEC.PY !! BRING HERE.
        P3D = (p3d * dk).sum() 
        sigmaRho = np.sqrt(P3D)

        
    if 0:
        alphaFile = open("sflaser_records.txt",'a')
        #alphaFile.write("\nThe R value %f"%(the_R))
        #alphaFile.write("\nThe R^(1/2) value %f"%(the_R_sqrt))
        alphaFile.write("\nSigma 3D from P3D %f"%(sigmaRho))
        #alphaFile.write("\nSigma 3D_2D %f"%(sigmaRho))
        alphaFile.close()

    return 
