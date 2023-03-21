
'''
calculates sigma rho straight from the data 
'''
from startlaser import * 


def Calculate(first_field,region,cg_field,cg_fieldtwo):
    
    # if we were not using a covergrid
    if first_field == 'density': 
        field_i = region['density']
        cell_volume = region['cell_volume']

        field_mean_i = np.mean(field_i)
        field_ave_i = ( (field_i*cell_volume).sum()/cell_volume.sum() )

    if 1:  # from the covergrid
        field_ii= cg_field.flatten()
        cv = cg_fieldtwo.flatten()

        field_mean_ii = np.mean(field_ii)
        field_ave_pii = ( (field_ii*cv).sum()/cv.sum() )
        # drop the units to be safe
        field_ave_ii = field_ave_pii.v 

    sfield_i = field_i/field_mean_i 
    sigma_sq_i = (((sfield_i - 1) ** 2)*cell_volume).sum() / cell_volume.sum()
    sigma_i = np.sqrt(sigma_sq_i)

    sfield_ii = field_ii/field_mean_ii 
    sigma_sq_ii = (((sfield_ii - 1) ** 2)*cv).sum() / cv.sum()
    # this could have also been coded as: 
    #sigma_sq_ii = (((field_ii - field_mean_ii) ** 2)*cv).sum() / cv.sum()

    sigma_ii = np.sqrt(sigma_sq_ii)

    if 0:
        pdb.set_trace()
    
    if 0:  # write to a file if you wish
        alphaFile = open("sflaser_records.txt",'a')
        alphaFile.write("\nsigma_rho calculated cg: %f"%sigma_ii)
        alphaFile.close()

    return 
