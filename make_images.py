
'''
    notes: 
    how does it look
'''
import matplotlib
import matplotlib.pyplot as plt
import pdb


def Make_Images(cg_projection,pj_axis):

    cg_proj = cg_projection.v
    #pdb.set_trace()
    
    plt.imshow(cg_proj,interpolation='nearest',origin='lower')
    plt.savefig('imgRho_proj0375_%d.png'%pj_axis)
    plt.close()

    return
