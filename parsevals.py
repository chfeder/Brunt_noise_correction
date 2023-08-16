
from startlaser import *

def parsevals_th(dim, p2d, p3d, dk, kvalues, sigmaN_sq=None, noisy=False, psn2d=0, psn3d=0):

    if noisy == True and dim == 1: 
        p2d = psn2d
    if noisy == True and dim == 2: 
        p3d = psn3d

    if dim == 1:
        if 1:
            P2D = ((p2d[1:]*dk).sum()) *dk
            P3D = ((p3d[1:]*dk).sum()) *dk
            sigmaRho_sq = sigmaN_sq * (P3D/P2D)
            sigmaRho = np.sqrt(sigmaRho_sq)

            the_R = P2D/P3D  
            the_R_sqrt = np.sqrt(P2D/P3D)

            if 0:  #print it out 
                #print('sigma ^2 from Brunts', sigmaRho_sq)
                print('sigma from Brunts', sigmaRho) 
                #print('the_R', the_R)                        
                print('the_R^sqrt', the_R_sqrt) 

        # testing the turbulence_driving_analysis script,def get_Brunt_sqrt_R(power_spectrum_map, k, exclude_k0=True) 
        if 0:
            kvalues.resize((217,217),refcheck=False)
            kvalues[1]=kvalues[0]
            p2dd = p2d.copy()  #otherwise it throws: ValueError: cannot resize this array: it does not own its data
            p2dd.resize((217,217),refcheck=False)
            p2d[1]=p2d[0]
            CF_sqrt_R = turbulence_driving_analysis.get_Brunt_sqrt_R(p2dd,kvalues)  
            sigmaRho_cf = np.sqrt(sigmaN_sq)/CF_sqrt_R   

            print('sigma cf',sigmaRho_cf) 
            print('the_R^sqrt cf',CF_sqrt_R)  
         

    if dim == 2: 
        P3D = ((p3d[1:] * dk).sum()) *dk   #the same if [1:] is taken out
        sigmaRho = np.sqrt(P3D)
        print('sigma from P3D', sigmaRho)

        
    if 0:  # write out what you wish
        alphaFile = open("sflaser_records.txt",'a')
        alphaFile.write("\nSigma 3D from P3D %f"%(sigmaRho))
        alphaFile.close()

    return sigmaRho 

