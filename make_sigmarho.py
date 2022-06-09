
'''
notes:
Fit to a function
Measure
'''
def Make_SigmaRho(bins,cbins,pdf):
    # Fit a gaussian or a gaussian+powerlaw
    # Sigma_log_rho = ( PDF* (log(rho)-mu)^2) .sum()
    #   Make sure this math is right.
    
    dx = bins[1:]-bins[:-1]
    norm = (pdf * dx).sum()
    meanX = (cbins * pdf * dx).sum()
    SigmaSquared = ((cbins - meanX)**2 * pdf * dx).sum()

    return
