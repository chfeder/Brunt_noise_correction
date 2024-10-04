
import platform
ver=platform.python_version()
python_version=int(ver[0])
if python_version == 3:
    print('python version 3...')
    from importlib import reload

import yt
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import ipdb
import math
import scipy
from scipy.optimize import curve_fit
from scipy import fft
import scipy.stats as stats
from numpy import inf 
import turbulence_driving_analysis
reload(turbulence_driving_analysis)

import cfpack as cfp
from cfpack import stop,print 
from cfpack import matplotlibrc   #should be okay to delete the matplotlibs above
import h5py


# the py file imports 
from yt_names import *
import extract_s1
reload(extract_s1)
import project
reload(project)
import parsevals
reload(parsevals)
import calculate
reload(calculate)
import noiseimage
reload(noiseimage)
#import make_images
#reload(make_images)


''' 
# to import things from sub directories
path_list= ["./the_sub_directory_1", "./the_subdirectory_2",...\
            ...]
for directory in path_list:
    if directory not in sys.path:
        sys.path += [directory]
'''
