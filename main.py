
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

'''
THE DATA: the simulations are in /scratch/ek9/ccf100/nif
* 9 different simulations with a few snapshots:
* NIF_hdf5_plt_cnt_0*
* frames = [000,125,250,275,300,325,350,375,425]
'''
the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_50mu_bubbles/NIF_hdf5_plt_cnt_0375')


'''
to obtain SIGMA_RHO:
if field_main = any velocity, then field_two = 'density', field_three = pressure
'''
field_main = 'density'
field_two = None 
field_three = None 
projection = True 
dimension = 1
proj_axis = 1  #or simply axis 
polytropic_index = 1 # +/-0.1, for 50um voids

the_level = 0
delta = 0.0175 
if 0:
    # trimmed analysis region, a box
    the_ledge = the_ds.arr([-delta,0.33,-delta],'code_length')
    the_redge = the_ds.arr([delta,0.365,delta],'code_length')
if 1:
    # original analysis region
    the_ledge = the_ds.arr([-0.06,0.33,-0.06],'code_length')
    the_redge = the_ds.arr([0.06,0.365,0.06],'code_length')

the_center =0.5*the_ledge + 0.5*the_redge
mindx = the_ds.index.get_smallest_dx()
size = ((the_redge - the_ledge)/mindx).astype('int') 
the_region = the_ds.region(the_center, the_ledge, the_redge)

section = extract_s1.extract_section_1(the_ds,the_level,the_ledge,size)
section_blur = motion_blur.Motion_Blur()  #tbd; simply returns
cg_first,cg_second,cg_third = project.Project(section,proj_axis,projection,field_main,field_two,field_three)
if 0:   
    images = make_images.Make_Images(cg_first,proj_axis) 
if 0: 
    calculation = calculate.Calculate(field_main,the_region,cg_first,cg_second) 
if 1: # this could be split into two definitions/files PDF and FIT
    bins, cbins, pdf, rhorhoavg = make_pdfNfit.Make_PDFnFIT(field_main,projection,proj_axis,cg_first,cg_second,cg_third) 
if 1:
    sigma_sq = make_sigmasquared.Make_SigmaSquared(field_main,bins,cbins,pdf,rhorhoavg)
if 1:
    powerspec2D, powerspec3D, infxyz = powerspec.PowerSpec(cg_first,dimension)
if 1:
    sigma_3D = parsevals.parsevals_th(dimension,sigma_sq, powerspec2D, powerspec3D, infxyz)


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
sigma_v = make_sigmavel.Make_SigmaVel(field_main,bins_v,cbins_v,pdf_v) 
'''


'''
TO OBTAIN b: 
'''
if 0: 
    obtain_b.Obtain_b()  #sigma_rho,sigma_v,polytropic_index)

