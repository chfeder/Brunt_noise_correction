
from startlaser import * 
'''
notes:
Fit to a function
Measure
'''

def Make_SigmaSquared(the_field,bins,cbins,pdf):
    
    ds = bins[1:]-bins[:-1]
    norm = (pdf * ds).sum()
    meanX = (cbins * pdf * ds).sum()


    if the_field == 'density':
        sigmasquared = ((np.exp(cbins) - 1)**2 * pdf * ds).sum()
        sigma = np.sqrt(sigmasquared)
    if the_field == 'overdensity':
        sigmasquared = (np.exp(cbins)**2 * pdf * ds).sum()  #since I logged the plotarray
        sigma = np.sqrt(sigmasquared)

    if 0:
        alphaFile = open("sflaser_records.txt",'a')
        #alphaFile.write("\nField %s, \nNorm %f, \nMeanX %f, \nSigma^2 %f"%(the_field,norm,meanX,sigmasquared))
        alphaFile.write("\nEq. 5 Sigma^2 of N/N_o %f"%(sigmasquared))
        #alphaFile.write("\nEq. 5 Sigma rho/rho_o %f"%(sigma))
        alphaFile.close()

    return sigmasquared
