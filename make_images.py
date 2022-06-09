
'''
    notes: 
    how does it look
'''
from main import *

def Make_Images(array):

    plt.imshow(array,interpolation='nearest',origin='lower')
    plt.savefig('imgRho_snap0325.png')
    plt.close()

    return
