
from startlaser import *
'''
NOTES:
1. Many layers of complexity.
Just project. cg.sum(axis=2)
Include optical effects
    Resolution of CCD
    PSF of the source
    Response of the CCD, Optical Depth effects
2. Opacity
Transmission = e^( -a Column )
3. Image = Source * Transmission
Add beam width
Add CCD resolution
Add CCD response


For later:
Proj = Convolve with Gaussian (need the PSF framing camera)
'''

def Project(covergrid,pj_axis,project,field_one,field_two=None,field_three=None):

        if field_one == 'overdensity': 
            rho_0 = (covergrid['density'] * covergrid['cell_volume']).sum()/covergrid['cell_volume'].sum()
            covergrid.set_field_parameter('rho_0',rho_0)

        velocities = ['sound_speed','velocity_x', 'velocity_y', 'velocity_z']

        if project == True:
            cg_main = covergrid[field_one].sum(axis=pj_axis)
            cg_two = None
            cg_three = None

        if project == False:
            cg_main = covergrid[field_one] 
            if field_two == 'cell_volume':
                cg_two = covergrid[field_two]
            else:
                cg_two = None
            cg_three = None

            '''
            if field_one in velocities[0]:  #think of a better way
                cg_main = covergrid['density'] 
                cg_two = covergrid[field_two]
                cg_three = covergrid[field_three]    
            if field_one in velocities[1:]:
                cg_main = covergrid[field_one]
                cg_two = covergrid[field_two]
                cg_three = covergrid[field_three]
            '''

        return cg_main, cg_two, cg_three
