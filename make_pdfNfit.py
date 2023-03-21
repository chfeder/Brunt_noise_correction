
'''
    notes: 
    PDF = histogram(3DField)
    and fit to a function

    this helps calculate sigma rho from the PDF
'''
from startlaser import *


def Make_PDFnFIT(the_field,project,axis,array_main,array_two=None,array_three=None):

    the_array = array_main.flatten() 
    the_meanarray = np.mean(the_array)  
    if array_two != None:
        the_sarray = array_two.flatten()
        rho = the_sarray 
    if array_three != None:
        the_tarray = array_three.flatten()
        p_th = the_tarray

    if the_field == 'density':
        the_rhorhoavg = the_array/the_meanarray
        the_plotarray = np.log(the_rhorhoavg)
    if the_field == 'overdensity':
        #pdb.set_trace()
        the_plotarray = np.log(abs(the_array))  #edited 


    velocities = ['sound_speed','velocity_x', 'velocity_y', 'velocity_z']
    if the_field in velocities: 
        gamma = 5/3
        c_s = np.sqrt((gamma * p_th)/rho)
        the_plotarray = (the_array - the_meanarray)/c_s 
    if the_field in velocities[0]:
        the_plotarray = c_s/np.mean(c_s)


    the_min = the_plotarray.min()
    the_max = the_plotarray.max()

    the_prebins = math.sqrt(len(the_plotarray))  #for isqrt python 3.8 is required
    the_bins = np.linspace(the_min, the_max,num=320)  #number of bins approximated to match paper

    the_PDF, xbins, patches = plt.hist(the_plotarray, bins=the_bins, histtype='step', weights = None, density=True,color='b')
    bin_centers = 0.5*(xbins[1:]+xbins[:-1])


    plot_x = bin_centers
    plot_y = the_PDF
    # THE CURVE
    if 0:
        plt.plot(plot_x,plot_y, c='k')

    # THE CURVE-FITS
    if the_field == 'density':
        def func_GaussPL(s,st,N,sigma,mu,p_o,a):
            y = np.zeros_like(plot_x)
            mask = s<st  #need a way to include st in the PL
            y[mask] = ( N/np.sqrt(2*np.pi*sigma**2) ) * np.exp(-( ((s[mask]-mu)**2)/(2*sigma**2) ) )
            y[~mask] = N*p_o*np.exp(-a*s[~mask])
            return y 
         
        if project == True:
            popt, pcov = curve_fit(func_GaussPL, xdata=plot_x, ydata=plot_y, p0 = [0.13,1,1,1,1,-5])
            print('popt',popt)
            plt.plot(plot_x, func_GaussPL(plot_x, *popt), 'r--')  #how to mark portions in diff colors
        if project == False:
            popt, pcov = curve_fit(func_GaussPL, xdata=plot_x, ydata=plot_y, p0 = [0.09,1,1,1,1,5])
            print('popt',popt)
            plt.plot(plot_x, func_GaussPL(plot_x, *popt), 'r--')  #how to mark portions in diff colors
        plt.ylabel('PDF')

    xyz = ['x','y','z']
    if the_field in velocities[1:]: 
        def func_gauss(xyz,sigma,mu): 
            y = ( 1/np.sqrt(2*np.pi*sigma**2) ) * np.exp(-( ((xyz-mu)**2)/(2*sigma**2) ) ) 
            return y
        popt, pcov = curve_fit(func_gauss, xdata=plot_x, ydata=plot_y)
        plt.plot(plot_x, func_gauss(plot_x, *popt), 'g--') 
        plt.xlabel(r'$\mathcal{M}_%s$'%xyz[axis])  #automate the axis labels...
        plt.ylabel('PDF')

    if the_field in velocities[0]:
        def func_PLGauss(s,st,A,B,sigma,mu,a):
            y = np.zeros_like(plot_x)
            mask = s<=st
            y[mask] = A*np.exp(a*s[mask])
            y[~mask] = ( B/np.sqrt(2*np.pi*sigma**2) ) * np.exp(-( ((s[~mask]-mu)**2)/(2*sigma**2) ) )
            return y 
        popt, pcov = curve_fit(func_PLGauss, xdata=plot_x, ydata=plot_y,p0=[0.95,1,1,1,1,-3])
        print('popt',popt)
        plt.plot(plot_x, func_PLGauss(plot_x, *popt), 'r--') 
        plt.xlabel(r'$c_s/ \left\langle c_s \right\rangle$')
        plt.ylabel('PDF')


    # POLISH
    plt.yscale('log')
    plt.ylim(1e-3,1e1)
    if the_field == 'density':
        plt.xlim(-0.5,0.8)
    if the_field in velocities[1:]: 
        plt.xlim(-1.5,1.5)
    if the_field == 'sound_speed':
        plt.xlim(0.2,1.5)

    plt.savefig('sort_plots/pdf%s_375_corrected.png'%(the_field))
    plt.close('all')
    return xbins, plot_x, plot_y 
          

if 0:  #DEBUG
    pdb.set_trace() 
