
'''
notes:
Fit to a function
Measure
'''
from main import *

def Make_SigmaRho(bins,cbins,pdf):
    # fit to a function: Fit a gaussian or a gaussian+powerlaw
    # measure: Sigma_log_rho = ( PDF* (log(rho)-mu)^2) .sum()
   

    dx = bins[1:]-bins[:-1]
    norm = (pdf * dx).sum()
    meanX = (cbins * pdf * dx).sum()
    sigmasquared = ((cbins - meanX)**2 * pdf * dx).sum()

    if 1:
        alphaFile = open("sflaser_records.txt",'a')
        alphaFile.write("Norm %f, \nMeanX %f, \nSigma^2 %f"%(norm,meanX,sigmasquared))
        alphaFile.close()



    return
