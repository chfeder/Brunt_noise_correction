
'''
notes:
Fit to a function
Measure
'''

def Make_SigmaSquared(the_field,bins,cbins,pdf):
    
    dx = bins[1:]-bins[:-1]
    norm = (pdf * dx).sum()
    meanX = (cbins * pdf * dx).sum()
    sigmasquared = ((cbins - meanX)**2 * pdf * dx).sum()

    if 0:
        alphaFile = open("sflaser_records.txt",'a')
        alphaFile.write("\nField %s, \nNorm %f, \nMeanX %f, \nSigma^2 %f"%(the_field,norm,meanX,sigmasquared))
        alphaFile.close()

    return sigmasquared
