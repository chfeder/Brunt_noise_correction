
'''
main.py: in progress 
'''
import yt
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import platform
ver=platform.python_version()
python_version=int(ver[0])
if python_version == 3:
    print('python version 3...')
    from importlib import reload
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
import make_sigmarho
reload(make_sigmarho)
import make_sigmavel
reload(make_sigmavel)
import produce_b
reload(produce_b)

'''
THE DATA: the simulations are in /scratch/ek9/ccf100/nif
* 9 different simulations with a few snapshots:
* NIF_hdf5_plt_cnt_0*
* frames = [000,125,250,275,300,325,350,375,425]
'''
the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_50mu_bubbles/NIF_hdf5_plt_cnt_0375')


'''
to obtain SIGMA_RHO:
'''
field_main = 'velocity_z'
field_two ='density' 
field_three ='pressure' 
proj_axis = 0  #doesnt make a difference if NOT making projection; assign it a None 

the_level = 0
the_ledge = the_ds.arr([-0.06,0.3,-0.06],'code_length')
the_redge = the_ds.arr([0.06,0.334,0.06],'code_length')
the_center =0.5*the_ledge + 0.5*the_redge
mindx = the_ds.index.get_smallest_dx()
size = ((the_redge - the_ledge)/mindx).astype('int') 


section = extract_s1.extract_section_1(the_ds,the_level,the_ledge,size)
section_blur = motion_blur.Motion_Blur()  #tbd; simply returns
cg_first,cg_second,cg_third = project.Project(section,proj_axis,field_main,field_two,field_three)
if 0:   
    images = make_images.Make_Images(cg_second)  #buggie 
if 0: # this could be split into two definitions/files PDF and FIT
bins, cbins, pdf = make_pdfNfit.Make_PDFnFIT(field_main,cg_first,cg_second,cg_third) 
if 0:
    sigma_rho = make_sigmarho.Make_SigmaRho(field_main,bins,cbins,pdf)


'''
to obtain SIGMA_V:
'''
'''
field = 'velocity'
proj_axis_vel = 2
#section_v = extract_s2.extract_section_2(the_ds2, the_level,the_ledge2,size2)  #for the more suitable region; using section for now
project_vel = project.Project(section,field,proj_axis_vel,the_center,the_redge,the_ledge)
#ProbablyOtherPhysicsTBD
bins_v, cbins_v, pdf_v = make_pdf.Make_PDF(project_vel,field) 
SigmaV = make_sigmavel.Make_SigmaVel(field_main,bins_v,cbins_v,pdf_v) 
'''
'''
TO OBTAIN b: 
'''
#dp_b = produce_b.Produce_B(SigmaRho,SigmaV)