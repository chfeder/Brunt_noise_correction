
'''
realize what moment maps are...
study: turbulence_driving_analysis.py
do: own idea & searches
'''
from startlaser import *

#import cf
#import cfplot as cfp  #...OHHH
#import cfpack as cfp  #pip install cfpack
# OR
# https://spectral-cube.readthedocs.io/en/latest/
from spectral_cube import SpectralCube

# from script learning
def moment_map(cg_main):

    cfp.plot_map(cg_main, log=True, cmap_label="column density", show=True)







