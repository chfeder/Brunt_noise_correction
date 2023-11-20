
from startlaser import * 
from collections import defaultdict

# for slices we need this encantation
from matplotlib import rcParams
rcParams['text.latex.preamble'] = r'\newcommand{\mathdefault}[1][]{}'

void_size = '100um'
if void_size=='12p5um': 
    the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_12p5mu_bubbles/NIF_hdf5_plt_cnt_0375')
if void_size=='50um': 
    the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_50mu_bubbles/NIF_hdf5_plt_cnt_0375')
if void_size=='100um':
    the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_100mu_bubbles/NIF_hdf5_plt_cnt_0375')

'''
to obtain SIGMA_RHO:
'''
# SETUP:
field_main = 'density'
field_two = 'cell_volume'
field_three = None 

polytropic_index = 1 # +/-0.1, for 50um voids
the_level = 0

# CONSTRUCTING THE VOLUMES, pick x or z:
box = 'sx_box'  

# M-sized boxes.
delta_m = 0.02  #(1.2/3)/2 
edge_mtop = 0.3675 
edge_mbot = 0.3275
viz_top = 0.4
viz_bot = 0.0

# variable volumes
delta = 0.01

# SIGNAL TO NOISE RATIO and iterationss: 
target = [inf]
many = 1


for pull in range(9): 
    pull = pull + 1
    if box == 'sx_box':
        the_ledge = the_ds.arr([-delta*pull, edge_mbot, -delta_m],'code_length') 
        the_redge = the_ds.arr([delta*pull, edge_mtop, delta_m],'code_length')  
        #the_ledge = the_ds.arr([-delta*pull, viz_bot, -delta_m],'code_length')  
        #he_redge = the_ds.arr([delta*pull, viz_top, delta_m],'code_length')  
    if box == 'sz_box':
        the_ledge = the_ds.arr([-delta_m, edge_mbot, -delta*pull],'code_length')
        the_redge = the_ds.arr([delta_m, edge_mtop, delta*pull],'code_length')
        #the_ledge = the_ds.arr([-delta_m, viz_bot, -delta*pull],'code_length')
        #the_redge = the_ds.arr([delta_m, viz_top, delta*pull],'code_length')
    
    the_center =0.5*the_ledge + 0.5*the_redge
    mindx = the_ds.index.get_smallest_dx()
    size = ((the_redge - the_ledge)/mindx).astype('int') 
    section = extract_s1.extract_section_1(the_ds,the_level,the_ledge,size)

    if pull==2:
        projection = True 
        dimension = 1
        proj_axis = 0  

        cg_first,_,_ = project.Project(section,proj_axis,projection,field_main,field_two,field_three) 
        _, p2d, p32d, dk32, ks = powerspec.powerspec(cg_first, dimension, target, many)    
        cfp.plot(p2d, ks, color='blue') #,ylog=True,ylim=[10e1,10e11], save='powerspec3d2d_%s_%s'%(pull,box))
        cfp.plot(p32d, ks, color='green') #,ylog=True,ylim=[10e1,10e11], save='powerspec3d2d_%s_%s'%(pull,box))

        if 1:
            projection = False 
            dimension = 2
            cg_first_i,_,_ = project.Project(section,proj_axis,projection,field_main,field_two,field_three) 
            _, _, p3d, dk3, ks = powerspec.powerspec(cg_first_i, dimension, target, many)    
            cfp.plot(p3d, ks, color='orange', ylog=True,ylim=[10e1,10e10], save='powerspec3d_3d2d_%s_%s_%s'%(pull,box,void_size))
            
        # same plots, not same values, focus is on the plots right now...
    