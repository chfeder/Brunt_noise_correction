
'''
    notes: 
    PDF = histogram( Projection )
'''
from main import *

def Make_PDF(frb):

    # take a peak
    #plt.imshow(np.log10(np.array(frb['gas','density'])),interpolation='nearest',origin='lower')
    #plt.savefig('imgRho_snap0325.png')

    the_x = np.log10(np.array(frb['gas','density']))
    the_weight = np.log10(np.array(frb['gas','cell_volume']))
    the_array, xbins = np.histogram(the_x, weights = None, density=True)
    bin_centers = 0.5*(xbins[1:]+xbins[:-1])
    plot_x = bin_centers
    plot_y = the_array
    plt.plot(plot_x,plot_y, c='k')

    plt.savefig('pdfRho_snap325.png')

    return pdf  #edit
          
