
'''
    notes: 
    how does it look
'''
from startlaser import *


def Make_Images(cg_projection,pj_axis):

    cg_proj = cg_projection.v
    
    plt.imshow(cg_proj,interpolation='nearest',origin='lower')
    plt.savefig('imgRho_proj0375_%d.png'%pj_axis)
    plt.close()

    return
