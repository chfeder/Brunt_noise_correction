
'''
main.py: in progress 
tbe: to be edited
'''
import yt
from extract_s1.py import *
from project.py import *
from make_pdf.py import *

'''
THE DATA: the simulations are in /scratch/ek9/ccf100/nif
* 9 different simulations with a few snapshots:
* NIF_hdf5_plt_cnt_0*
* frames = [000,125,250,275,300,32,350,375,425]
'''
the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_50mu_bubbles/NIF_hdf5_plt_cnt_0325')


'''
to obtain SIGMA_RHO
'''
the_level = 0
the_leftedge = [-1.5,-2.0,-1.5]

section = extract_section_1( the_ds,the_level,the_leftedge )
#section1 = Motion_Blur( section )  #tbe
projection = Project( section )  #tbe
pdf = Make_PDF( projection )
SigmaRho = MakeSigmaRho( PDF )


'''
to obtain SIGMA_V
'''
Section_v = extract_section_2( argument )
ProbablyOtherPhysicsTBD
SigmaV = MakeSigmaV( Section_v )


'''
TO OBTAIN b 
'''
B = ProduceB( SigmaRho, SigmaV )
