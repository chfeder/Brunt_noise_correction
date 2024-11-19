
'''
NOTES:
1. Many layers of complexity.
Just project. cg.sum(axis=2)
Include optical effects
    Resolution of CCD
    PSF of the source
    Response of the CCD, Optical Depth effects
2. Opacity
Transmission = e^( -a Column )
3. Image = Source * Transmission
Add beam width
Add CCD resolution
Add CCD response


For later:
Proj = Convolve with Gaussian (need the PSF framing camera)
'''
from startlaser import *

def Project(covergrid, pj_axis, project, field_one, field_two=None):

    if project == True:
        pre_cg_main_sum = covergrid[field_one].sum(axis=pj_axis)

        cg_main_mean = pre_cg_main_sum.mean()
        cg_main_norm = pre_cg_main_sum/cg_main_mean

        cg_main = pre_cg_main_sum
        
        cg_two = covergrid[field_two].sum(axis=pj_axis)

    if project == False:
        cg_main = covergrid[field_one] 
        if field_two == 'cell_volume':
            cg_two = covergrid[field_two]
        else:
            cg_two = None

    return cg_main, cg_main_norm, cg_two
