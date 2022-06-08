
'''
main.py: in progress 
tbe: to be edited
'''
import yt
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import platform
ver=platform.python_version()
python_version=int(ver[0])
if python_version == 3:
    from importlib import reload
import extract_s1
reload(extract_s1)
import motion_blur
reload(motion_blur)
import project
reload(project)
import make_pdf
reload(make_pdf)

'''
THE DATA: the simulations are in /scratch/ek9/ccf100/nif
* 9 different simulations with a few snapshots:
* NIF_hdf5_plt_cnt_0*
* frames = [000,125,250,275,300,32,350,375,425]
'''
the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_50mu_bubbles/NIF_hdf5_plt_cnt_0325')


'''
to obtain SIGMA_RHO:
'''
the_level = 0
the_ledge = the_ds.arr([-0.06,0.3,-0.06],'code_length')
the_redge = the_ds.arr([0.06,0.334,0.06],'code_length')
the_center =0.5*the_ledge + 0.5*the_redge
proj_axis = 1
mindx = the_ds.index.get_smallest_dx()
size = ((the_redge - the_ledge)/mindx).astype('int') 


section = extract_s1.extract_section_1(the_ds,the_level,the_ledge,size)
section1 = motion_blur.Motion_Blur()  #tbe; simply returns
projection = project.Project(section,proj_axis,the_center,the_redge,the_ledge)  #revise
pdf = make_pdf.Make_PDF(projection)  #revise  
#SigmaRho = MakeSigmaRho(pdf)


'''
to obtain SIGMA_V:

Section_v = extract_section_2( argument )
ProbablyOtherPhysicsTBD
SigmaV = MakeSigmaV( Section_v )
'''

'''
TO OBTAIN b: 

B = ProduceB( SigmaRho, SigmaV )
'''

