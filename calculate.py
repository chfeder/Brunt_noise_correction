
'''
calculates sigma rho straight from the data 
'''
from startlaser import * 


def Calculate(first_field, cg_field=None, cg_fieldtwo=None):
    
    if 1:  # from the covergrid
        field_ii= cg_field.flatten()
        cv = cg_fieldtwo.flatten()

        field_mean_ii = np.mean(field_ii)
        #field_ave_pii = ( (field_ii*cv).sum()/cv.sum() )
        field_ave_ii = ( (field_ii*cv).sum()/cv.sum() )
        # drop the units to be safe
        #field_ave_ii = field_ave_pii.v 

    sfield_ii = field_ii/field_mean_ii 
    sigma_sq_ii = (((sfield_ii - 1) ** 2)*cv).sum() / cv.sum()
    if 0:
        print('sigma_2_calc', sigma_sq_ii)
    # this could have also been coded as: 
    #sigma_sq_ii = (((field_ii - field_mean_ii) ** 2)*cv).sum() / cv.sum()
    sigma_ii = np.sqrt(sigma_sq_ii)

    if 0: 
        print('sigma 2D calculated', sigma_ii) 
    if 0:  # write to a file if you wish
        alphaFile = open("sflaser_records.txt",'a')
        alphaFile.write("\nsigma_rho calculated cg: %f"%sigma_ii)
        alphaFile.close()

    return sigma_sq_ii, sigma_ii 
