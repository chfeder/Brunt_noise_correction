

# the typical imports, e.g. starter1.py
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
import pdb
import math
import scipy
from scipy.optimize import curve_fit
from scipy import fft
import scipy.stats as stats


# the not so typical imports, e.g. starter2.py
from yt_names import *
import extract_s1
reload(extract_s1)
import motion_blur
reload(motion_blur)
import project
reload(project)
import make_images
reload(make_images)
import make_pdfNfit
reload(make_pdfNfit)
import make_sigmasquared
reload(make_sigmasquared)
import make_sigmavel
reload(make_sigmavel)
import parsevals
reload(parsevals)
import obtain_b
reload(obtain_b)
import calculate
reload(calculate)
import powerspec
reload(powerspec)

# in parsevals.py only...
#import turbulence_driving_analysis
#reload(turbulence_driving_analysis)


''' 
# to import things from sub directories
path_list= ["./the_sub_directory_1", "./the_subdirectory_2",...\
            ...]
for directory in path_list:
    if directory not in sys.path:
        sys.path += [directory]
'''
