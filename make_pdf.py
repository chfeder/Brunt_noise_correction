
'''
    notes: 
    PDF = histogram( Projection )
'''
from main import *

def Make_PDF(array):

    the_PDF, xbins = np.histogram(array, weights = None, density=True)
    bin_centers = 0.5*(xbins[1:]+xbins[:-1])
    plot_x = bin_centers
    plot_y = the_PDF
    plt.plot(plot_x,plot_y, c='k')

    plt.savefig('pdfRho_snap325.png')

    return xbins,plot_x,plot_y 
          
