
import numpy as np
'''
notes:
'''

def Obtain_b():#sigmasqrho,sigmasqvel,p_index):
    # compute the driving parameter
    # the arguments passed in squared values

    # TEMPORARY FROM SFLASER RECORDS:
    p_index = 0.2
    sigmasqrho = 0.053404
    sigmasqvel = 0.048430 + 0.085634 + 0.047294

    sigma_rho = np.sqrt(sigmasqrho)
    mach_num = sigmasqvel

    # EQUATION 1
    driving_param_one = np.sqrt(sigmasqrho/sigmasqvel)
    # EQUATION 10
    numerator = (sigma_rho**(2*p_index)) - 1
    denominator = mach_num * p_index * (1-(1/(sigma_rho**2))) 
    driving_param_ten = np.sqrt(numerator/denominator) 
    print(driving_param_one,driving_param_ten)

    if 1:
        alphaFile = open("sflaser_records.txt",'a')
        #alphaFile.write("\ndrive param eq 1: %f"%driving_param_one)
        alphaFile.write("\ndrive param eq 10: %f"%driving_param_ten)
        alphaFile.close()

    return
