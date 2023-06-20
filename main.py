
from startlaser import * 
'''
main.py: in progress 

THE DATA: the simulations are in /scratch/ek9/ccf100/nif
* 9 different simulations with a few snapshots:
* NIF_hdf5_plt_cnt_0*
* frames = [000,125,250,275,300,325,350,375,425]
'''

the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_50mu_bubbles/NIF_hdf5_plt_cnt_0375')

'''
to obtain SIGMA_RHO:
if field_main = any velocity, then field_two = 'density', field_three = pressure
derived field (an option): to take the mean out, aka k=0 mode
'''

def add_field(thing):
    def overdensity(field,data):
        rho_0 = data.get_field_parameter('rho_0')
        #rho_ratio = (data['density']/rho_0).v
        #return (rho_ratio-1)
        # OR SAME THING:
        return(((data['density'] - rho_0)/rho_0))
    thing.add_field(YT_overdensity,function=overdensity)#,units='dimensionless')
add_field(the_ds)


field_main = 'density'
field_two = 'cell_volume'
field_three = None 
projection = True 
dimension =1
proj_axis = 1  #or simply axis 
polytropic_index = 1 # +/-0.1, for 50um voids
the_level = 0

# which box? right, left, center, not a box
box = 'c_box'
# target signal to noise suggested ratios: inf,100(20db),10(10db), 5(6.9897db), 3(4.7712db), 1(0db)
target_snr_db = [100, 20, 10, 6.9897, 4.7712, 0, -10, -20, -30, -40, -50, -60, -70, -80, -90, -100]
noise = True
noisy_powerspec2D = 0
noisy_powerspec3D = 0
target_snr = 0  #since np.sqrt() accepts value 'only' up to 19 ooms!

''' 
quick next: triple check code
NEXT! pseudo-
for left, center, right box positions:
    for i in target_snr:
        get sigma rho 3D from 2D
        plot sigma3D_2D vs snr_db
'''

delta = 0.0175 
if box == 'c_box':
    # trimmed analysis region, a box
    the_ledge = the_ds.arr([-delta,0.33,-delta],'code_length')
    the_redge = the_ds.arr([delta,0.365,delta],'code_length')
if box == 'r_box':
    the_ledge = the_ds.arr([delta,0.33,delta],'code_length')
    the_redge = the_ds.arr([delta+2*delta,0.365,delta+2*delta],'code_length')
if box == 'l_box':
    the_ledge = the_ds.arr([-delta-2*delta,0.33,-delta-2*delta],'code_length')
    the_redge = the_ds.arr([-delta,0.365,-delta],'code_length')
if box == 'notabox':
    # original analysis region
    the_ledge = the_ds.arr([-0.06,0.33,-0.06],'code_length')
    the_redge = the_ds.arr([0.06,0.365,0.06],'code_length')

the_center =0.5*the_ledge + 0.5*the_redge
mindx = the_ds.index.get_smallest_dx()
size = ((the_redge - the_ledge)/mindx).astype('int') 
the_region = the_ds.region(the_center,the_ledge, the_redge)

section = extract_s1.extract_section_1(the_ds,the_level,the_ledge,size)
section_blur = motion_blur.Motion_Blur()  #tbd; simply returns
cg_first,cg_second,cg_third = project.Project(section,proj_axis,projection,field_main,field_two,field_three)

if 0:   
    images = make_images.Make_Images(cg_first,proj_axis) 
if 0: 
    calculate.Calculate(field_main,the_region,cg_first,cg_second) 

if 1:
    # the following could be split into two definitions/files PDF and FIT
    bins, cbins, pdf = make_pdfNfit.Make_PDFnFIT(field_main,projection,proj_axis,cg_first,cg_second,cg_third) 
    sigma_sq = make_sigmasquared.Make_SigmaSquared(field_main,bins,cbins,pdf)
    if projection == False: 
        sigma_sq = None
if 1:
    four_amps, powerspec2D, powerspec3D, infxyz, ks = powerspec.PowerSpec(cg_first,cg_second,dimension)
if 1:
    noisy_powerspec2D, noisy_powerspec3D = signaltonoise.signaltonoise(cg_first, dimension,target_snr,four_amps,powerspec2D,powerspec3D) 
if 1:
    sigma_3D = parsevals.parsevals_th(dimension,powerspec2D,powerspec3D,infxyz,ks,sigma_sq,noise,noisy_powerspec2D,noisy_powerspec3D) 


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

