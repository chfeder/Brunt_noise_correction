
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
'''
# For later:
# Proj = Convolve with Gaussian (need the PSF framing camera)

def Project(covergrid,pj_axis,project,field_one,field_two=None,field_three=None):

        velocities = ['sound_speed','velocity_x', 'velocity_y', 'velocity_z']

        if project == True:
            cg_main = covergrid[field_one].sum(axis=pj_axis)
            cg_two = None
            cg_three = None

        if project == False:
            if field_one == 'density':
                cg_main = covergrid[field_one]
                cg_two = None
                cg_three = None
            if field_one in velocities[0]:  #think of a better way
                cg_main = covergrid['density'] 
                cg_two = covergrid[field_two]
                cg_three = covergrid[field_three]    
            if field_one in velocities[1:]:
                cg_main = covergrid[field_one]
                cg_two = covergrid[field_two]
                cg_three = covergrid[field_three]

        return cg_main, cg_two, cg_three

#let;s see if the passing of field_one/two works...