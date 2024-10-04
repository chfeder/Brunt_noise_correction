
'''
calculates sigma rho straight from the data 
'''
from startlaser import * 


def Calculate(first_field, cg_field=None, cg_fieldtwo=None):
    
    field= cg_field.flatten()
    cv = cg_fieldtwo.flatten()

    field_mean = np.mean(field)

    sfield = field/field_mean 
    sigma_sq_ii = (((sfield - 1) ** 2)*cv).sum() / cv.sum()
    sigma_ii = np.sqrt(sigma_sq_ii)

    return sigma_sq_ii, sigma_ii 


    # NOTE:
    # this could have also been coded as: 
    #sigma_sq_ii = (((field_ii - field_mean_ii) ** 2)*cv).sum() / cv.sum()