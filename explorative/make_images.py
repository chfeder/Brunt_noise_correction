
'''
    notes: 
    how does it look
'''
from startlaser import *


def Make_Images(ds, cg_projection, pj_axis, box, center, width):

    #cg_proj = cg_projection.v
    
    if 0:
        plt.imshow(cg_proj,interpolation='nearest',origin='lower')
        plt.savefig('imgRho_proj0375_%d_%s.png'%(pj_axis,box))
        plt.close()

    if 1:
        projz = yt.SlicePlot(ds, pj_axis, ('gas','density'), data_source = cg_projection) #center=center, width=width)
        outname = 'density_slice0375_i_%d_%s.png'%(pj_axis,box)
        projz.save()
        print('saved ',box)

    return
