
'''
    notes: 
    how does it look
'''

def Make_Images(array):

    plt.imshow(array,interpolation='nearest',origin='lower')
    plt.savefig('imgRho_snap0375.png')
    plt.close()

    return
