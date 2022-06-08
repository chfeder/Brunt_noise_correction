
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
'''
# For later:
# Proj = Convolve with Gausian (need the PSF framing camera)
from main import *

def Project(covergrid,pj_axis,center,redge,ledge):
    Proj = covergrid['gas','density'].sum(axis=pj_axis)

    return Proj
