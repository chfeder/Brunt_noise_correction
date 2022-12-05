
import numpy as np
'''
notes:
'''

def Obtain_b():#sigmasqrho,sigmasqvel,p_index):
    # compute the driving parameter
    # the arguments passed in squared values

    # TEMPORARY FROM SFLASER RECORDS:
    p_index = 1
    sigmasqrho = 0.084967
    sigmasqvel = 0.025556 + 0.063517 + 0.024683

    sigma_rho = np.sqrt(sigmasqrho)
    mach_num = sigmasqvel

    numerator = (sigma_rho**(2*p_index)) - 1
    denominator = mach_num * p_index * (1-(1/(sigma_rho**2))) 
    driving_param = np.sqrt(numerator/denominator) 
    print(driving_param)

    if 1:
        alphaFile = open("sflaser_records.txt",'a')
        alphaFile.write("\ndrive param: %f"%driving_param)
        alphaFile.close()

    return
