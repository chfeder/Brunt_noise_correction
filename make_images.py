
import matplotlib
import matplotlib.pyplot as plt
'''
    notes: 
    how does it look
'''

def Make_Images(array,pj_axis):

    plt.imshow(array,interpolation='nearest',origin='lower')
    plt.savefig('imgRho_snap0375_%d.png'%pj_axis)
    plt.close()

    return
