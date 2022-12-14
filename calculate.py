

'''
calculate sigma rho from data as opposed to from PDFs
not equation (4), even more standard

'''
import matplotlib.pyplot as plt
import numpy as np
import pdb


def Calculate(first_field):

    field = first_field.flatten()
    field_ave = np.mean(field)
    sigma_sq = ((field - field_ave) ** 2).sum()
    sigma = np.sqrt(sigma_sq)

    # hm, but it should be of rho/rho_avg...

    if 0:
        pdb.set_trace()
    
    if 1:
        alphaFile = open("sflaser_records.txt",'a')
        alphaFile.write("\nsigma calculated: %f"%sigma)
        alphaFile.close()

    return
