
import platform
ver=platform.python_version()
python_version=int(ver[0])
if python_version == 3:
    print('python version 3...')
    from importlib import reload

import yt
from yt.funcs import mylog
mylog.setLevel(30)
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
import h5py

import turbulence_driving_analysis
reload(turbulence_driving_analysis)

import cfpack as cfp
from cfpack import stop, print
cfp.load_plot_style()

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
import fake_powerlaw
reload(fake_powerlaw)
