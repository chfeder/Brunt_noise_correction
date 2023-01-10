
'''
notes:
Fit to a function
Measure
'''
import pdb
import numpy as np

def Make_SigmaSquared(the_field,bins,cbins,pdf,prefield):
    
    dx = bins[1:]-bins[:-1]
    norm = (pdf * dx).sum()
    meanX = (cbins * pdf * dx).sum()
    if 0:
        sigmasquared = ((cbins - meanX)**2 * pdf * dx).sum()

    if 1:  #equation 5. 
        if the_field == 'density':
            sigmasquared_c = ((np.exp(cbins) - 1)**2 * pdf * dx).sum()
        else:
            sigmasquared_c = 0
        sigma = np.sqrt(sigmasquared_c)

    if 0:
        alphaFile = open("sflaser_records.txt",'a')
        #alphaFile.write("\nField %s, \nNorm %f, \nMeanX %f, \nSigma^2 %f"%(the_field,norm,meanX,sigmasquared))
        #alphaFile.write("\nEq. 5 Sigma^2 of N/N_o %f"%(sigmasquared_c))
        alphaFile.write("\nEq. 5 Sigma rho/rho_o %f"%(sigma))
        alphaFile.close()

    return sigmasquared_c
