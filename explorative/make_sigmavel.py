
'''
Section is a little later in time, a little further down the tube
SigmaV = Section[‘velocity_y’].std()
We’ll do something smarter when we grok get visar a little better
'''


def Make_SigmaVel(the_field,bins,cbins,pdf):
   
    dx = bins[1:]-bins[:-1]
    norm = (pdf * dx).sum()
    meanX = (cbins * pdf * dx).sum()
    sigmasquared = ((cbins - meanX)**2 * pdf * dx).sum()

    if 1:
        alphaFile = open("sflaser_records.txt",'a')
        alphaFile.write("\n %s \nNorm %f, \nMeanX %f, \nSigma^2 %f"%(the_field,norm,meanX,sigmasquared))
        alphaFile.close()

    return sigmasquared
