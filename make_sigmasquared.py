
'''
notes:
Fit to a function
Measure
'''

def Make_SigmaSquared(the_field,bins,cbins,pdf,prefield):
    
    dx = bins[1:]-bins[:-1]
    norm = (pdf * dx).sum()
    meanX = (cbins * pdf * dx).sum()
    sigmasquared = ((cbins - meanX)**2 * pdf * dx).sum()
 
    if the_field == 'density':
        sigmasquared_c = ((prefield - 1)**2 * pdf * dx).sum()
    else:
        sigmasquared_c = 0

    if 1:
        alphaFile = open("sflaser_records.txt",'a')
        alphaFile.write("\nField %s, \nNorm %f, \nMeanX %f, \nSigma^2 %f"%(the_field,norm,meanX,sigmasquared))
        alphaFile.write("\nCalculated Sigma^2 %f"%(sigmasquared_c))
        alphaFile.close()

    return sigmasquared
