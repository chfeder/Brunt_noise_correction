import yt
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import unyt
from unyt import cm, s


# 9 different simulations with a few snapshots:
# NIF_hdf5_plt_cnt_0*
#frames = [000,125,250,275,300,32,350,375,425]
ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_50mu_bubbles/NIF_hdf5_plt_cnt_0325')

# FOR NEXT PROJS, ETC, DO QUICKER WITH
#for fr in frames:
#    ds_list = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_50mu_bubbles/NIF_hdf5_plt_cnt_%04d'%frames[fr])
'''
HM: ds.region takes code length, explore...
SEE PAPER, also:
print(ds.domain_left_edge)
print(ds.domain_right_edge)
print(ds.domain_right_edge - ds.domain_left_edge)
'''
left_edge = ds.arr([-0.06,0.3,-0.06],'code_length') 
right_edge = ds.arr([0.06,0.334,0.06],'code_length')
center =0.5*left_edge + 0.5*right_edge
the_region = ds.region(center, left_edge, right_edge)#, fields=None, ds=None, field_parameters=None, data_source=None)

# check out the region
if 1:
    proj_axis = 2
    proj_d = ds.proj(('gas','density'),proj_axis,data_source=the_region)
    pw = proj_d.to_pw()
    pw.save()

    ix = ds.coordinates.x_axis[proj_axis]
    iy = ds.coordinates.y_axis[proj_axis]

    midpt = center 
    width = right_edge - left_edge
    width_2D = width[ix],width[iy]  
    #dx_inv = ds.domain_dimensions/ds.domain_width  #are the zones cubed?  #dir(ds), 
    #num_zones = width * dx_inv
    
    nx = np.unique(the_region['x']).size
    ny = np.unique(the_region['y']).size
    nz = np.unique(the_region['z']).size
    num_zones = nx, ny, nz
    rez = [num_zones[ix],num_zones[iy]] #[768,1024]
    
    frb_d = proj_d.to_frb(width_2D,rez,center=midpt)
    frb = frb_d
    plt.imshow(np.log10(np.array(frb['gas','density'])),interpolation='nearest',origin='lower')
    plt.savefig('imgRho_snap0325.png')
    raise

if 0:
    the_x = np.log10(np.array(frb['gas','density']))
    the_weight = np.log10(np.array(frb['gas','cell_volume']))
    the_array, xbins = np.histogram(the_x, weights = None, density=True)
    bin_centers = 0.5*(xbins[1:]+xbins[:-1])
    plot_x = bin_centers
    plot_y = the_array
    plt.plot(plot_x,plot_y, c='k')
    plt.savefig('pdfRho_snap325.png')

# TAKE 1 Calculation
if 0:
    # sigma_rho/rho_avg = b mach_number; mach_number = sigma_v/c_s
    # from Dhawalikar et al Nov2021, for 50um voids, <c_s> = 20km/s
    gamma = 5/3
    p_th = the_region['pressure']
    rho = the_region['density']
    rho_avg = np.mean(rho)
    c_s = np.sqrt((gamma * p_th)/rho)
    cs_avg = np.mean(c_s)
    #cs_avg = [2e6] * cm/s  #could have also done ds.array([#],'units')
  
    vx = the_region['velocity_x']
    vx_avg = np.mean(vx)  # np.mean() is ok as long as region remains constant
    vy = the_region['velocity_y']
    vy_avg = np.mean(vy)
    vz = the_region['velocity_z']
    vz_avg = np.mean(vz)
    v2 = vx**2 + vy**2 + vz**2
    v = np.sqrt(v2)

    sigma_v = np.std(v)
    sigma_rho = np.std(rho)

    # OR ..
    sigma_vel2 = (vx-vx_avg)**2 + (vy-vy_avg)**2 + (vz-vz_avg)**2 
    sigma_vel = np.sqrt(sigma_vel2)

    LHS = sigma_rho/rho_avg
    RHS = sigma_v/cs_avg
    b = LHS/RHS
    print('b ',b)  #dump numbers in a text file later

# quick check
if 0:
    projz = yt.ProjectionPlot(ds,'z',('gas','density'))
    projz.save()
    projy = yt.ProjectionPlot(ds,'y',('gas','density'))
    projy.save()
    projx = yt.ProjectionPlot(ds,'x',('gas','density'))
    projx.save()

