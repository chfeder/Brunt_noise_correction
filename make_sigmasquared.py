
'''
notes:
this calculates sigma rho from the PDF
'''
from startlaser import * 


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

    if 0:  # write out anything you wish
        alphaFile = open("sflaser_records.txt",'a')
        alphaFile.write("\nEq. 5 Sigma^2 of N/N_o %f"%(sigmasquared))
        alphaFile.close()

    return sigmasquared
