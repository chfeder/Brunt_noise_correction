
from startlaser import * 
from collections import defaultdict
import matplotlib.ticker as ticker
import matplotlib.lines as mlines
# GET INFO FROM STORAGE 
import pdb 


# POSSIBLE FIT FUNCTIONS
# discovery fit functions
def linfunct(x, a, b):  
    y = a * x + b 
    return y
def lininvfunct_og(x, a, b=0.25):
    y = a * x**(b) - 0.2 
    return y
def discovery(x, a, b, c):
    # originally an exponential function
    pre_y = 1 + x**-a 
    y = np.log(pre_y)
    return y
def tanfunct_simple(x, first, second): 
    y = first * np.tanh(second*x) + (1-first)  
    return y
def gaussianlike(x, a, b): 
    #y = a * np.exp(-(x-1)**2 / (2 * b**2)) + c
    y = a * np.exp(-(x)**2 /b) + 1
    return y
def sigmarho_fit_old_old(x,a):  
    pre_y = (lininvfunct(x, a, b=1) * tanfunct_simple(x, first=first, second=second)) 
    y = np.log(pre_y)
    return y 
def sigmarho_fit_old(x, a, b, c):
    y = 1/(a * x**(b)+c)
    return y 

# the fitting functions we actually used
def csnr_1(x, a, b):
    y = np.log( (x**a + b) / x**a )
    return y 
def csnr_2(x, c, d):
    y = np.log( x**c / (x**c + d) )
    return y 

def csnr(x, a, b, c, d):  
    cfactor = csnr_1(x,a,b) + csnr_2(x,c,d)
    return cfactor 


# GET THE FITTING 
def fitthebox(whichfigure, eldata, secondparam=0):  
    '''
    GET the median, the 16th ('-' correction), and the 84th ('+' correction) of fit params a, b, c, d, etc
    this then allows us to report param +/- correction 
    '''
    param_a = []
    param_b = []
    param_c = []
    param_d = [] 
    for j in range(numofboxes):  
        xdata = [0.1, 1, 3, 5, 10, 100] 
        oneboxproj = []
        for k in range(len(xdata)):
            prexdata = eldata[k][j]
            oneboxproj = np.append(oneboxproj, prexdata) 
        
        mask = ~np.isnan(oneboxproj)
        filt_xdata = np.array(xdata)[mask]
        filt_oneboxproj = oneboxproj[mask]  
        if whichfigure == 'sigmasdenoised':
            secondparam ='justtwo'
            print('inside the sigmasdenoised')
            params = {"a": [0.1, 1, 3.0], "b": [0.1, 1, 3.0]}
            #arr_safe = np.nan_to_num(oneboxproj, nan=1)
            rets = cfp.fit(csnr_1, filt_xdata, np.log(filt_oneboxproj), params=params)
        if whichfigure == 'sigmas':
            print('inside the sigmas')
            secondparam ='justtwo'
            params = {"c": [0.1, 1, 3.0], "d": [0.1, 1, 3.0]}
            arr_safe = np.nan_to_num(oneboxproj, nan=1)
            rets = cfp.fit(csnr_2, filt_xdata, np.log(filt_oneboxproj), params=params)
        if whichfigure == 'sigmastwo':
            print('inside the sigmas')
            secondparam ='gauss'
            params = {"a": [0.1, 1, 3.0], "b": [0.1, 1, 3.0], "c":[0.1, 1, 3.0],"d": [0.1, 1, 3.0] } 
            arr_safe = np.nan_to_num(oneboxproj, nan=1)
            rets = cfp.fit(csnr, filt_xdata, np.log(filt_oneboxproj), params=params)

        param_a = np.append(param_a, rets.popt[0])
        if secondparam == 'gauss' or secondparam == 'justtwo':
            param_b = np.append(param_b, rets.popt[1])
        if secondparam== 'gauss':
            param_c = np.append(param_c, rets.popt[2])
            param_d = np.append(param_d, rets.popt[3])

    med_a = np.percentile(param_a, 50) 
    p16_a = np.percentile(param_a, 16)
    p84_a = np.percentile(param_a, 84)
    if secondparam == 'gauss' or secondparam =='justtwo':
        med_b = np.percentile(param_b, 50)
        p16_b = np.percentile(param_b, 16)
        p84_b = np.percentile(param_b, 84)
    if secondparam == 'gauss':
        med_c = np.percentile(param_c, 50)
        med_d = np.percentile(param_d, 50)
        p16_c = np.percentile(param_c, 16)
        p16_d = np.percentile(param_d, 16)
        p84_c = np.percentile(param_c, 84)
        p84_d = np.percentile(param_d, 84)

    if secondparam ==0:
        med_b = p16_b = p84_b = 0
    if secondparam ==0 or secondparam == 'justtwo': 
        med_c = p16_c = p84_c = 0
        med_d = p16_d = p84_d = 0

    if secondparam == 0:
        param_meds = [med_a]
        param_16s = [p16_a] 
        param_84s = [p84_a] 
    if secondparam == 'justtwo':
        param_meds = [med_a, med_b]
        param_16s = [p16_a, p16_b]
        param_84s = [p84_a, p84_b]
    if secondparam == 'gauss':
        param_meds = [med_a, med_b, med_c, med_d]
        param_16s = [p16_a, p16_b, p16_c, med_d]
        param_84s = [p84_a, p84_b, p84_c, med_d]
    
    return param_meds, param_16s, param_84s


def plotthefitting(the_xxx, thefigures, themedians, sixteenths=None, eightyfours=None, fitting='no'):
    """ plot the fitting """
    if thefigures == 'sigmasdenoised':
        the_yrecipemed = csnr_1(the_xxx, *themedians)  
        if fitting=='yes': 
            the_yrecipelow = csnr_1(the_xxx, *sixteenths)
            the_yrecipeupp = csnr_1(the_xxx, *eightyfours)
    if thefigures == 'sigmas':
        the_yrecipemed = csnr_2(the_xxx, *themedians)  
        if fitting=='yes': 
            the_yrecipelow = csnr_2(the_xxx, *sixteenths)
            the_yrecipeupp = csnr_2(the_xxx, *eightyfours)
    if thefigures == 'sigmastwo':
        the_yrecipemed = csnr(the_xxx, *themedians)  
        if fitting=='yes': 
            the_yrecipelow = csnr(the_xxx, *sixteenths)
            the_yrecipeupp = csnr(the_xxx, *eightyfours)
    if fitting =='no':
        the_yrecipelow = 0  
        the_yrecipeupp = 0 
    return the_yrecipemed, the_yrecipelow, the_yrecipeupp

if 'kludge' not in dir():
    kludge ={}
def errorbarsplot(thefigure, ladata, void_name=None, secondfig=None, seconddata=None, void_errname=None, styleplot='singular'):
    """ plot the data """
    thefigs = [thefigure, secondfig]
    thedata = [ladata, seconddata]
    colores = [colors[draw], colors[draw], colors[draw], colors[draw]]
    labels = [void_name, void_errname]  #this used to be singular no "s"
    #alphas = [0.7, 0.3, 0.4, 0.5]
    alphas = [0.6, 0.6, 0.6, 0.6]
    markers = ['s', 'o', '*', 'P']
    for num, figs in enumerate(thefigs):
        the_medians = []
        the_uppers = []
        the_lowers = []
        for k in range(len(thedata[num])):
            if secondfig:
                for m in range(len(thedata[num][k])):
                    if seconddata[k][m] == 0:
                        thedata[num][k][m] = np.nan
            the_medians.append(np.nanpercentile(thedata[num][k], 50))
            the_lowers.append(np.nanpercentile(thedata[num][k], 16))
            the_uppers.append(np.nanpercentile(thedata[num][k], 84))
        the_medians = np.array(the_medians)
        the_lowers = np.array(the_lowers)
        the_uppers = np.array(the_uppers)
        lower_error = the_medians - the_lowers
        upper_error = the_uppers - the_medians
        yerrs = [lower_error[1:], upper_error[1:]] 

        if styleplot =='singular': 
            if typeofsigma == 'free' and thefigure=='sigmas':
                ax.errorbar(xdata, the_medians, yerr=[lower_error, upper_error], linestyle=None, fmt=markers[1], \
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5, label=labels[num])
            else:
                ax.errorbar(xdata[1:], the_medians[1:], yerr=yerrs, linestyle=None, fmt=markers[draw], \
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5, label=labels[num])

        if styleplot =='subpanel':
            if thefigure =='sigmasdenoised' and typeofsigma=='free':
                ax1.errorbar(xdata[1:], the_medians[1:], yerr=yerrs, linestyle=None, fmt=markers[draw], \
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5, label=labels[num])
                if void_name =='12.5$\mathrm{\mu m}$':
                    ax1.axhline(y=the_medians[k], linestyle='dashed', color='gray', alpha=0.7, label=r'$\sigma_{\rho/\rho_0\mathrm{Brunt,denoised}}\mathrm{(SNR\!\to\!\infty)}$')
                else:
                    ax1.axhline(y=the_medians[k], linestyle='dashed', color='gray', alpha=0.7)

            if thefigure =='sigmasdenoised' and typeofsigma==None:
                ax2.errorbar(xdata[1:], the_medians[1:], yerr=yerrs, linestyle=None, fmt=markers[draw], \
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5)
            if thefigure =='sigmas': 
                ax2.errorbar(xdata[1:], the_medians[1:], yerr=yerrs, linestyle=None, fmt=markers[draw], \
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5, label=labels[num])
            if thefigure =='sigmasdenoised' and typeofsigma=='cthree':
                if panels=='two':
                    ax2.errorbar(xdata[1:], the_medians[1:], yerr=yerrs, linestyle=None, fmt=markers[draw], \
                                ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5)
                if panels=='three':
                    ax3.errorbar(xdata[1:], the_medians[1:], yerr=yerrs, linestyle=None, fmt=markers[draw], \
                                ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5)

        if styleplot =='subpanel':
            if thefigure =='rsqrtfull':
                print('inside plotting rhalf plot ',void_name)
                ax1.errorbar(xdata[1:], the_medians[1:], yerr=yerrs, linestyle=None, fmt=markers[draw], \
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5, label=labels[num])
            if thefigure =='rsqrttrue':
                ax2.errorbar(xdata[1:], the_medians[1:], yerr=yerrs, linestyle=None, fmt=markers[draw], \
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5, label=labels[num])
                if void_name =='12.5$\mathrm{\mu m}$':
                    ax2.axhline(y=the_medians[k], linestyle='dashed', color='gray', alpha=0.7, label=r'$R^{1/2}_{\mathrm{denoised}}\mathrm{(SNR\!\to\!\infty)}$')
                else:
                    ax2.axhline(y=the_medians[k], linestyle='dashed', color='gray', alpha=0.7)
        if secondfig == None:
            break
        if secondfig:
            break

    # FITTING THE PLOTS
    # WARNING: add which figures can enter the fitting.
    if thefitting == 'one' or thefitting == 'two' or thefitting == 'three':
        if typeofsigma=='cthree':
            thefigure='sigmastwo'

        # to fit each box
        the_xx = np.linspace(xdata[1:].min(), xdata[1:].max(), num=2500)  

        # call the fittings procedure definition
        if thefigure != 'sigmastwo' or fit=='yes': 
            meds, lows, upps = fitthebox(thefigure, ladata[1:])
            if thefitting == 'one': 
                medianos1.append(meds)
            if thefitting == 'two': 
                medianos2.append(meds)
            if thefitting == 'three': 
                medianos3.append(meds)

        if void_name == '100$\mathrm{\mu m}$':
            if thefitting == 'one':
                mediano_1 = np.mean(medianos1,axis=0) 
                ln_ymed, ln_ylow, ln_yupp = plotthefitting(the_xx, thefigure, mediano_1)
                #ymed1 = np.exp(ln_ymed) 
                #kludge['ymed1'] = ymed1
            if thefitting == 'two':
                mediano_2 = np.mean(medianos2,axis=0) 
                ln_ymed, ln_ylow, ln_yupp = plotthefitting(the_xx, thefigure, mediano_2)
                #ymed2 = np.exp(ln_ymed) 
                #kludge['ymed2'] = ymed2
            if thefitting == 'three':
                if fit =='no':
                    mediano_1 = np.mean(medianos1,axis=0) 
                    mediano_2 = np.mean(medianos2,axis=0) 
                    mediano = np.append(mediano_1, mediano_2)
                    print('the four params ',mediano)
                    ln_ymed, ln_ylow, ln_yupp = plotthefitting(the_xx, thefigure, mediano)
                    #ymed3 = kludge['ymed1']*kludge['ymed2']  
                    #ymedreal = np.exp(ln_ymed) 
                if fit =='yes':
                    mediano = np.mean(medianos3,axis=0) 
                    ln_ymed, ln_ylow, ln_yupp = plotthefitting(the_xx, thefigure, mediano)
            ymed = np.exp(ln_ymed)

            if styleplot=='subpanel':
                if whichfigtwo !='sigmasdenoised' and typeofsigma==None:
                    ax2.plot(the_xx, ymed, color='grey', linestyle='dashed', linewidth=3, label='$C_1\mathrm{(SNR)}$ fit', alpha=0.7)
                    #ax2.plot(the_xx, kludge['ymed1'], color='grey', linestyle='dashed', linewidth=3, label='$C_1\mathrm{(SNR)}$ fit', alpha=0.7)
                if whichfigtwo=='sigmasdenoised':  
                    ax2.plot(the_xx, ymed, color='grey', linestyle='dashed', linewidth=3, label='$C_2\mathrm{(SNR)}$ fit', alpha=0.7)
                    #ax2.plot(the_xx, kludge['ymed2'], color='k', linestyle='dashed', linewidth=3, label='$C_2\mathrm{(SNR)}$ fit', alpha=0.5)
                if typeofsigma=='cthree':
                    if panels=='two':
                        ax2.plot(the_xx, ymed, color='r', linestyle='dashed', linewidth=3, label='$C\mathrm{(SNR)}$', alpha=0.7)
                    if panels=='three':
                        ax3.plot(the_xx, ymed, color='r', linestyle='--', linewidth=3, label='$C\mathrm{(SNR)}$', alpha=0.7)
                        #ax3.plot(the_xx, ymed3, color='r', linestyle='--', linewidth=3, label='$C\mathrm{(SNR)}$', alpha=0.7)
    return  


def plotfig(plotstyle):
    if plotstyle == 'singular':
        if whichfig == 'kvals':
            label = r'$k_{\mathrm{noise}}/k_{\mathrm{max}}$'
            ylims = (1.5e-2,1.5e0)
            yscales = 'log'
            ax.set_xlabel("SNR", fontsize=13)
            ax.set_ylabel(label, fontsize=13)
            ax.legend(loc='lower right', fontsize=12, labelspacing=0.2)
        if whichfig == 'rsqrtfull':
            label = r'$R^{1/2}$'
            ylims = None
            yscales = 'linear'
        if whichfig == 'rsqrttrue':
            label = r'$R^{1/2}_{\mathrm{denoised}}$'
            ylims = None
            yscales = 'linear'
        if whichfig == 'ratios':
            label = r'$R^{1/2}_{\mathrm{denoised}}/R^{1/2}$'
            ylims = None
            yscales = 'linear'
        if whichfig == 'sigmas' and typeofsigma==None:
            label = 'Noise correction C(SNR)'
            ax.set_xlabel("SNR", fontsize=12)
            ax.set_ylabel(label, fontsize=12)
            ylims = None
            yscales = 'log'
            handles, labels = ax.get_legend_handles_labels()
            desired_label = 'fit' 
            if desired_label in labels:
                index = labels.index(desired_label)
                desired_handle = handles.pop(index)
                desired_label = labels.pop(index)
                handles.append(desired_handle)
                labels.append(desired_label)
                ax.legend(handles, labels, loc='lower right', fontsize=12, labelspacing=0.2) 
        if typeofsigma == 'free' and whichfig =='sigmas':
            label = r'$\sigma_{\mathrm{\rho/\rho_o,Brunt}}$' #\, is a thin space, \: a medium space  
            ax.set_xlabel("SNR", fontsize=13)
            ax.set_ylabel(label, fontsize=13)
            ax.legend(loc='best', fontsize=12)
            ylims = (1e-1,1.5e2)
            yscales = 'log'
        if typeofsigma =='free' and whichfig =='sigmasdenoised':
            yscales = 'linear'
            ylims = (0,1)
            label = r'$\sigma_{\mathrm{\rho/\rho_o,Brunt,denoised}}$'
            ax.set_xlabel("SNR", fontsize=13)
            ax.set_ylabel(label, fontsize=13)
            ax.legend(loc='upper left', fontsize=12)
        if whichfig == 'sigmasdenoised' and typeofsigma==None:
            label = 'Noise correction C(SNR)'
            ax.set_xlabel("SNR", fontsize=13)
            ax.set_ylabel(label, fontsize=13)
            ylims = None#(0.5, 3)#(1e-1,0.5e1)
            yscales = 'log'
        if whichfig == 'sigmas' and typeofsigma=='cthree':
            label = 'Noise correction C(SNR)'
            ax.set_xlabel("SNR", fontsize=13)
            ax.set_ylabel(label, fontsize=13)
            ylims = None#(0.5, 3)#(1e-1,0.5e1)
            yscales = 'log'
        if whichfig == 'columnddenoise':
            label = r'$\sigma_{N/N_0,\mathrm{denoised}}$'
            ylims = (0,0.3)
            yscales = 'linear'
            ax.set_xlabel("SNR", fontsize=13)
            ax.set_ylabel(label, fontsize=13)

        #yscales='log'  #WARNING, TEMP FIX
        #ylims=None     #WARNING, TEMP FIX
        ax.set(xscale='log', yscale=yscales, ylim=ylims)

    if plotstyle == 'subpanel' and whichfig!= 'sigmasdenoised':
        ax1.set(xscale='log', yscale='linear') 
        y1label=r'$R^{1/2}$'
        ax1.set_ylabel(y1label, fontsize=21)
        ax1.tick_params(axis='both', which='major', labelsize=21)

        ax2.set(xscale='log', yscale='linear') 
        y2label=r'$R^{1/2}_{\mathrm{denoised}}$'
        ax2.set_ylabel(y2label, fontsize=21)
        ax2.tick_params(axis='both', which='major', labelsize=21)
        ax2.set_xlabel("SNR", fontsize=21)

        handles, labels = ax2.get_legend_handles_labels()
        desired_label = r'$R^{1/2}_{\mathrm{denoised}}\mathrm{(SNR\!\to\!\infty)}$'
        if desired_label in labels:
            index = labels.index(desired_label)
            desired_handle = handles.pop(index)
            desired_label = labels.pop(index)
            handles.append(desired_handle)
            labels.append(desired_label)
            ax2.legend(handles, labels, loc='best', fontsize=20, labelspacing=0.2) 


    if plotstyle == 'subpanel' and whichfig== 'sigmasdenoised':
        ax1.set(xscale='log', yscale='linear') 
        ylims=(0,1)
        ax1.set(ylim=ylims)
        y1label = r'$\sigma_{\mathrm{\rho/\rho_o,Brunt,denoised}}$'
        ax1.set_ylabel(y1label, fontsize=21)
        ax1.tick_params(axis='both', which='major', labelsize=21)

        handles, labels = ax1.get_legend_handles_labels()
        desired_label = r'$\sigma_{\rho/\rho_0\mathrm{Brunt,denoised}}\mathrm{(SNR\!\to\!\infty)}$'
        if desired_label in labels:
            index = labels.index(desired_label)
            desired_handle = handles.pop(index)
            desired_label = labels.pop(index)
            handles.append(desired_handle)
            labels.append(desired_label)
            ax1.legend(handles, labels, loc='upper left', fontsize=20, labelspacing=0.2) 

        ax2.set(xscale='log', yscale='log') 
        y2label = 'Noise correction $C_i\mathrm{(SNR)}$'
        ax2.set_ylabel(y2label, fontsize=21)
        ax2.tick_params(axis='both', which='major', labelsize=21)
        ax2.legend(loc='best', fontsize=21)
        if panels=='two':
            ax2.set_xlabel('SNR', fontsize=21)

        if panels=='three':
            ax3.set(xscale='log', yscale='log') 
            y3label = 'Noise correction $C\mathrm{(SNR)}$'
            ax3.set_ylabel(y3label, fontsize=21)
            ax3.tick_params(axis='both', which='major', labelsize=21)
            ax3.legend(loc='lower right', fontsize=21)
            ax3.set_xlabel('SNR', fontsize=21)


    # really saving now
    fig.tight_layout()   #trouble for sigmas?
    fig.savefig('sort_plots/cubes15_spectra15_%s_%s_i1knobox.pdf'%(whichfig,powerlaw))
    print('%s_saved!'%whichfig)
    plt.clf()




# MAIN: LOOP FOR THE THREE VOIDS!
# figure options: sigmas(typeofsigma='free',singular), kvals(singular,typeofsigma=None), rsqrtfull (subpanel, two),
#                 sigmasdenoised(subpanel, three); kvals2D(typeofsigma=None)
numofboxes = 15
whichfig = 'kvals'
ave = False  # there's normalized data by the average and non-normalized data
whichfigtwo = None
typeofsigma = None 
powerlaw='two'
target_snr = [0.01, 0.1, 1, 3, 5, 10, 100, inf]  
xdata = np.array(target_snr[:-1])

k_plot = 'divkmax'  # divkmax or knoise
plotstyle = 'singular'  # singular or subpanel
panels = 'three'  # two or three
fit = 'no'  #does not fit the final correction factor (as it should be)
thefitting=None

voids = ['12p5mu', '50mu', '100mu']
voids_formal = ['12.5$\mathrm{\mu m}$', '50$\mathrm{\mu m}$', '100$\mathrm{\mu m}$']
voids_eformal = ['$12.5_{\mu m, \mathrm{denoised}}$', '$50_{\mu m, \mathrm{denoised}}$', '$100_{\mu m, \mathrm{denoised}}$',]
colors = ['magenta', 'b', 'g']

if typeofsigma == 'free':  #or 'free' for the one void figure 
    voids = ['50mu']  
    voids_formal = ['50$\mathrm{\mu m}$']
    colors = ['b']    

if 1:
    if whichfig != 'kvals2D' and plotstyle=='singular':
        fig,ax = plt.subplots(1,1) 
    if whichfig == 'rsqrtfull' and plotstyle=='subpanel':
        fig,(ax1, ax2) = plt.subplots(2,1, sharex=True, figsize=(7, 9))  

    if whichfig == 'sigmasdenoised' and plotstyle=='subpanel':
        if panels=='two':
            fig,(ax1, ax2) = plt.subplots(2,1, sharex=True, figsize=(7, 9))  
        if panels=='three':
            fig,(ax1, ax2, ax3) = plt.subplots(3,1, sharex=True, figsize=(8, 12))  
        medianos1=[]
        medianos2=[]
        medianos3=[]

    for draw, void_size in enumerate(voids):
        print('void size is ',void_size)
        #hfivefile = '/home/552/lv9897/sflaser/h5files/boxplot15_powspec14_i1k_half_%s.h5'%void_size
        hfivefile = '/home/552/lv9897/sflaser/h5files/15cubes15spectra_i1k_two_%s.h5'%void_size
        Fptr = h5py.File(hfivefile,'r')

        # SCATTERS + FITS
        if whichfig != 'kvals2D':

            # SINGULAR PLOTS
            if plotstyle == 'singular':
                data = list(Fptr['%s_%s'%(whichfig,void_size)])  #returns a list 
                print('inside singular')

                if whichfig == 'kvals' and k_plot=='divkmax':
                    for val, kvalues in enumerate(data):
                        data[val] = data[val]/data[-1] 
                    errorbarsplot(whichfig, data, void_name=voids_formal[draw])  
                if whichfig == 'rsqrttrue':
                    errorbarsplot(whichfig, data, void_name=voids_formal[draw], void_errname=voids_eformal[draw])  
                if whichfig == 'rsqrtfull':
                    errorbarsplot(whichfig, data, void_name=voids_formal[draw], void_errname=voids_eformal[draw])  
                if whichfig == 'ratios':
                    print('inside fig')
                    theparams = 1  #two params
                    errorbarsplot(whichfig, data, void_name=voids_formal[draw])  
                if whichfig == 'sigmas' and whichfigtwo != 'sigmasdenoised':
                    if typeofsigma == None:
                        infdata = data[-1]/data 
                        errorbarsplot(whichfig, infdata, void_name=voids_formal[draw])  
                    if typeofsigma == 'free':
                        errorbarsplot(whichfig, data, void_name=voids_formal[draw])  
                # NEW.
                if whichfig == 'columnddenoise':
                    errorbarsplot(whichfig, data, void_name=voids_formal[draw])  

                if whichfig == 'sigmasdenoised':
                    if typeofsigma == 'free':
                        errorbarsplot(whichfig, data, void_name=voids_formal[draw])  
                    if typeofsigma == None:
                        infdata = data[-1]/data 
                        errorbarsplot(whichfig, infdata, void_name=voids_formal[draw])  

                if whichfig == 'sigmas' and whichfigtwo=='sigmasdenoised':
                    figtwodata =list(Fptr['%s_%s'%(whichfigtwo,void_size)])  
                    if typeofsigma =='cthree':
                        # WARNING, THIS DOES NOT WORK YET
                        csnrdata1 = figtwodata[-1]/figtwodata
                        whichfig = 'sigmasdenoised'
                        thefitting='one'
                        errorbarsplot(whichfig, csnrdata1, void_name=voids_formal[draw], styleplot=plotstyle)  
                        csnrdata2 = np.array(figtwodata)/np.array(data)
                        whichfig = 'sigmas'
                        thefitting='two'
                        errorbarsplot(whichfig, csnrdata2, void_name=voids_formal[draw], styleplot=plotstyle)  
                        csnrdata_three = figtwodata[-1]/data
                        whichfig = 'sigmasdenoised'
                        thefitting = 'three'
                        errorbarsplot(whichfig, csnrdata_three, void_name=voids_formal[draw], styleplot=plotstyle)  


            # MULTIPANEL PLOTS
            if plotstyle == 'subpanel' and whichfig=='sigmasdenoised':
                data = list(Fptr['%s_%s'%(whichfig,void_size)])  #returns a list 
                typeofsigma ='free'
                thefitting=None
                errorbarsplot(whichfig, data, void_name=voids_formal[draw], styleplot=plotstyle)  

                typeofsigma = None
                infdata = data[-1]/data 
                thefitting='one'
                errorbarsplot(whichfig, infdata, void_name=voids_formal[draw], styleplot=plotstyle)  

                whichfig = 'sigmas'
                whichfigtwo = 'sigmasdenoised'
                typeofsigma = None
                data_sigmas = list(Fptr['%s_%s'%(whichfig,void_size)])  
                csnrdata = np.array(data)/np.array(data_sigmas)
                thefitting='two'
                errorbarsplot(whichfig, csnrdata, void_name=voids_formal[draw], styleplot=plotstyle)  

                whichfig='sigmasdenoised'
                whichfigtwo = None
                figtwodata = None
                typeofsigma = 'cthree'
                csnrdata_three = data[-1]/data_sigmas
                thefitting='three'
                errorbarsplot(whichfig, csnrdata_three, void_name=voids_formal[draw], styleplot=plotstyle)
                # note: saving happens post powerspectra plots


            if plotstyle == 'subpanel' and whichfig=='rsqrtfull':
                thefitting=None
                data = list(Fptr['%s_%s'%(whichfig,void_size)])  #returns a list 
                errorbarsplot(whichfig, data, void_name=voids_formal[draw], void_errname=voids_eformal[draw], styleplot=plotstyle)
                print('after rhalf plot ',void_size)

                whichfig = 'rsqrttrue'
                data = list(Fptr['%s_%s'%(whichfig,void_size)])  #returns a list 
                errorbarsplot(whichfig, data, void_name=voids_formal[draw], void_errname=voids_eformal[draw], styleplot=plotstyle)
                #              secondfig=whichfigtwo, seconddata=figtwodata)
                print('after true plot ',void_size)
                whichfig = 'rsqrtfull'


        # POWERSPECTRA PLOTS 
        if whichfig == 'kvals2D':
            figfour = 'kvals'
            figfive = 'powintersects'
            fig,ax = plt.subplots(1,1)
            
            whichfigtwo = 'thepowspec' 
            whichfigthree = 'thenoisepowspec'  
            ps_xdata = Fptr['%s_%s'%(whichfig,void_size)][()]  #returns a numpy array
            if k_plot == 'divkmax':  #for k/k_max
                for val, kvalues in enumerate(ps_xdata):
                    ps_xdata[val] = ps_xdata[val]/ps_xdata[-1] 
            if whichfigtwo:
                ydata = Fptr['%s_%s'%(whichfigtwo,void_size)][()]  
            if whichfigthree:
                nydata = Fptr['%s_%s'%(whichfigthree,void_size)][()]  

            target_snr = [1000000, 0.01, 0.1, 1, 3, 5, 10, 100]  
            x_offset = 0.05 #0.08, 05
            y_offset = 5500 
            y_offset100 = 9500 

            # MEDIANS ONLY
            if 1:  
                medians_snr = []
                medians_noise = []
                the_lowers = []
                the_uppers = []
                for iii in range(len(ydata)):
                    medians_snr.append(np.median(ydata[iii], axis=0))
                    medians_noise.append(np.median(nydata[iii], axis=0))
                    the_lowers.append(np.percentile(ydata[iii], 16, axis=0))
                    the_uppers.append(np.percentile(ydata[iii], 84, axis=0))
                #stop()

                for i, a in enumerate(target_snr):  #len of 7
                    if i == 0:
                        ax.plot(ps_xdata, medians_snr[i], color='r', alpha=1.0, linewidth=0.8)
                        continue
                    if i != len(medians_snr):
                        ax.plot(ps_xdata[1:], medians_snr[i][1:], color=colors[draw], alpha=1.0, linewidth=0.8)
                        ax.plot(ps_xdata[1:], medians_noise[i][1:], color='grey',alpha=1.0, linewidth=0.8, linestyle='dashed')
                        if a == 3:
                            plt.text(ps_xdata[-1] + x_offset, medians_noise[i][-1]+0.006, str(a), fontsize=12, ha='left', va='center') 
                            continue
                        if a == 5:
                            if void_size == '100mu':
                                plt.text(ps_xdata[-1] + x_offset, medians_noise[i][-1]-0.002, str(a), fontsize=12, ha='left', va='center') 
                            else:
                                plt.text(ps_xdata[-1] + x_offset, medians_noise[i][-1], str(a), fontsize=12, ha='left', va='center') 
                            continue
                        if a == 10:
                            if void_size == '100mu':
                                plt.text(ps_xdata[-1] + x_offset, medians_noise[i][-1]-0.0004, str(a), fontsize=12, ha='left', va='center') 
                            else:
                                plt.text(ps_xdata[-1] + x_offset, medians_noise[i][-1]-0.0001, str(a), fontsize=12, ha='left', va='center') 
                        else:
                            plt.text(ps_xdata[-1] + x_offset, medians_noise[i][-1], str(a), fontsize=12, ha='left', va='center') 

                    if i == len(target_snr[1:]): 
                        ax.plot(ps_xdata[1:], medians_snr[i][1:], color=colors[draw], alpha=1.0, linewidth=0.8, label='S+N') 
                        ax.plot(ps_xdata[1:], medians_noise[i][1:], color='grey',alpha=1.0, linewidth=0.8, linestyle='dashed', label='N') 
                        ax.plot(ps_xdata[1:], medians_snr[0][1:], color='r', alpha=1.0, linewidth=0.8, label='S') 
                        ax.fill_between(ps_xdata[1:], the_lowers[0][1:], the_uppers[0][1:], color='r', alpha=0.2)
                    #if i == 1: 
                    #    plt.text(ps_xdata[-1]+x_offset, medians_noise[i][-1]+y_offset100, 'SNR', fontsize=12, ha='left', va='center')  

                data = list(Fptr['%s_%s'%(figfour,void_size)])  #returns a list 
                data_p = list(Fptr['%s_%s'%(figfive,void_size)]) 
                for val, kvalues in enumerate(data):
                    data[val] = data[val]/data[-1] 
                yrange = np.linspace(0, 1, 10000)
                for i, an in enumerate(target_snr[2:7]):  #len of 7
                   i = i+1
                   the_medians = np.nanpercentile(data[i], 50)
                   the_medians_low = np.nanpercentile(data[i], 16)
                   the_medians_high = np.nanpercentile(data[i], 84)
                   the_medians_ints = np.nanpercentile(data_p[i], 50) 
                   the_medians_intshigh = np.nanpercentile(data_p[i], 84) 

                   ax.plot([the_medians,the_medians], [0,the_medians_ints], color=colors[draw], linestyle='dotted', alpha=0.8)
                   ax.fill_betweenx(yrange, the_medians_low, the_medians_high, color=colors[draw], alpha=0.1, \
                                    where=(yrange<=the_medians_intshigh), linestyle='dotted')

            # SAVING! 
            if k_plot == 'divkmax':  #for k/k_max
                labeling = r'$k/k_{\mathrm{max}}$'
            if ave==False:
                ax.set(xscale='log', yscale= 'log', xlim=(2e-2,1.5e0), ylim= (1e-7,3.5e4+y_offset))  #(1e-2,1.7e0+x_offset) 
            if ave==True:
                ax.set(xscale='log', yscale= 'log', xlim=(2e-2,1.5e0), ylim= None)  #(1e-2,1.7e0+x_offset) 
            ylabel=r'$P_{\mathrm{2D}}$'
            ax.set_xlabel(labeling, fontsize=12)
            ax.set_ylabel(ylabel, fontsize=12)

            if void_size == '100mu': 
                main_legend = ax.legend(loc='upper left', ncol=1, fontsize=12, labelspacing=0.1, bbox_to_anchor=(0.01, 1.02))
            else:
                main_legend = ax.legend(loc='upper left', ncol=1, fontsize=12, labelspacing=0.1, bbox_to_anchor=(0.01, 1.01))
            ax.add_artist(main_legend)
            custom_line = mlines.Line2D([], [], color=colors[draw], linestyle='dotted', alpha= 0.8, label=r'$k_{\mathrm{noise}}/k_{\mathrm{max}}$')
            k_legend = ax.legend(handles=[custom_line], loc='upper left', fontsize=12, bbox_to_anchor=(0.23,1.01))

            plt.savefig('sort_plots/cubes15_spectra15_i1k_%s_%s_%s.pdf'%(powerlaw,k_plot,void_size))
            print('saved spectra14!')
            plt.clf()
            Fptr.close()         

            # A QUICK CHECK TO SEE IF P2D VS K LOG-LOG POWERLAW IS INDEED 1. YES. 
            if 0:
                quickx = np.log10(ps_xdata)
                quicky = np.log10(nydata)
                rets = cfp.fit(linfunct, quickx[1:], quicky[-1][1:])
                stop()  #fix.
        Fptr.close() 


    if whichfig != 'kvals2D':
        plotfig(plotstyle) 


# Figure 7 with noise spec k=1, k=1/2, k=0, k=2, one void size only
# WARNING: changed lines 162 - 166
if 0:
    void_size = '12p5mu'
    #voids_formal = ['50$\mathrm{\mu m}$, pl=1', '50$\mathrm{\mu m}$, pl=1/2']
    voids_formal = ['pl=1', 'pl=1/2', 'pl=0', 'pl=2']
    plotstyle='subpanel'
    whichfig = 'sigmasdenoised' 
    hfivefile_one = '/home/552/lv9897/sflaser/h5files/15cubes15spectra_i1k_one_%s.h5'%void_size
    hfivefile_half = '/home/552/lv9897/sflaser/h5files/15cubes15spectra_i1k_half_%s.h5'%void_size
    hfivefile_zero = '/home/552/lv9897/sflaser/h5files/15cubes15spectra_i1k_zero_%s.h5'%void_size
    hfivefile_two = '/home/552/lv9897/sflaser/h5files/15cubes15spectra_i1k_two_%s.h5'%void_size
    hfivefiles = [hfivefile_one, hfivefile_half, hfivefile_zero, hfivefile_two]    
    fig,(ax1, ax2, ax3) = plt.subplots(3,1, sharex=True, figsize=(8, 12))  
    medianos1=[]
    medianos2=[]
    medianos3=[]
    colors=['b','cyan','g','orange']
    for draw, hfivefile in enumerate(hfivefiles):    
        Fptr = h5py.File(hfivefile,'r') 
        data = list(Fptr['%s_%s'%(whichfig,void_size)])  #returns a list 
        typeofsigma ='free'
        thefitting=None
        errorbarsplot(whichfig, data, void_name=voids_formal[draw], styleplot=plotstyle)  

        print('draw',draw)

        typeofsigma = None
        infdata = data[-1]/data 
        thefitting='one'
        errorbarsplot(whichfig, infdata, void_name=voids_formal[draw], styleplot=plotstyle)  

        whichfig = 'sigmas'
        whichfigtwo = 'sigmasdenoised'
        typeofsigma = None
        data_sigmas = list(Fptr['%s_%s'%(whichfig,void_size)])  
        csnrdata = np.array(data)/np.array(data_sigmas)
        thefitting='two'
        errorbarsplot(whichfig, csnrdata, void_name=voids_formal[draw], styleplot=plotstyle)  

        whichfig='sigmasdenoised'
        whichfigtwo = None
        figtwodata = None
        typeofsigma = 'cthree'
        csnrdata_three = data[-1]/data_sigmas
        thefitting='three'
        errorbarsplot(whichfig, csnrdata_three, void_name=voids_formal[draw], styleplot=plotstyle)
    plotfig(plotstyle)
    
# plotting sigma rho calculated vs sigma rho brunt 
if 0:
    boxes = ['l_box', 'c_box', 'r_box', 'b_box', 'f_box']
    voids = ['50mu', '100mu', '12p5mu']
    colores = ['b', 'g', 'orange']
    valuesfoam = []
    def percent_difference_3(value1, value2, value3, value4, value5):
        mean_value = (value1 + value2 + value3 + value4 + value5) / 5 
        pd1 = abs(value1 - mean_value) / mean_value * 100
        pd2 = abs(value2 - mean_value) / mean_value * 100
        pd3 = abs(value3 - mean_value) / mean_value * 100 
        pd4 = abs(value4 - mean_value) / mean_value * 100 
        pd5 = abs(value5 - mean_value) / mean_value * 100 
        return pd1, pd2, pd3, pd4, pd5
    for val, boxesin in enumerate(boxes):
        values = []
        for num, void_size in enumerate(voids):
            hfivefile = '/home/552/lv9897/sflaser/h5files/15sigma3D_%s.h5'%void_size
            hfivefile_B = '/home/552/lv9897/sflaser/h5files/15nonoise_%s.h5'%void_size
            Fptr = h5py.File(hfivefile,'r')
            Fptr_B = h5py.File(hfivefile_B,'r')
            ydata_3D = list(Fptr['sigma_3Ds_%s'%void_size])  
            ydata_Brunt = list(Fptr_B['sigma_brunts_%s'%void_size])  
            ydata_2Ds = list(Fptr_B['sigma_2Ds_%s'%void_size])  
            ydata_rhalf = list(Fptr_B['sigma_rhalfs_%s'%void_size])  

            values.append(ydata_3D[val])
            if void_size == '12p5mu':
                valuesfoam.append(values)
                stop()
        if 0:
            porcentajes = []
            def relative_error(trueval, brunti):
                pdiff = (np.abs(brunti - trueval) / trueval) * 100
                return pdiff
            for i in range(len(ydata_3D)):
                pdiff = relative_error(ydata_3D[i], ydata_Brunt[i])
                porcentajes.append(pdiff)
            mean_pdiff = np.mean(np.array(porcentajes))
        xdata = np.linspace(1,15,15)
        fig,ax = plt.subplots(1,1)
        ax.scatter(xdata, ydata_Brunt, c=colores[num], alpha=0.7, label='brunt')
        ax.scatter(xdata, ydata_3D, c='k', alpha=0.5, label='3D')
        ax.set(ylim=(0.0,1.0), xlabel='number of boxes', ylabel=r'$\sigma_{\rho}$') 
        ax.legend(loc='best')
        print('saving, saved')
        fig.savefig('sort_plots/15boxes_brunt_3D_%s'%void_size)
        plt.close()
    Fptr.close() 



'''
# MAKE A BOXPLOT AND SCATTER PLOT - left for learning purposes; how to hack a boxplot. 
def boxplotscatter(thefigure, ladata, theparams=0, secondfig=None, seconddata=None):
    fig,ax = plt.subplots(1,1)
    # for the plot within the figure
    if thefigure in {'kvals', 'rsqrtfull'}:
        right, top, width, height = [0.6, 0.2, 0.25, 0.25]
        ax2 = fig.add_axes([right, top, width, height])
    if thefigure in {'ratios', 'sigmas'}:
        right, bottom, width, height = [0.6, 0.6, 0.25, 0.25]
        ax2 = fig.add_axes([right, bottom, width, height])

    # funny order to make the boxplots with a plot overlayed to work
    ax.set_xscale('log')
    ax.set_xlim(np.array(subxdata).min()/1.20, np.array(subxdata).max()*1.20)
    ax2.set_xscale('log')
    ax2.set_xlim(np.array(xdata).min()/1.20, np.array(xdata).max()*1.20)

    # the boxplot!  WARNING: this sets only for 'sigmas' option 
    if thefigure == 'sigmas':
        the_data = ladata/ladata[-1]
        data_lista = the_data.tolist()
    else:
        data_lista = ladata
    subboxplot = ax.boxplot(data_lista[2:], showfliers=True, showmeans=False, labels=subxdata, positions=subxdata) 
    boxplot = ax2.boxplot(data_lista, showfliers=True, showmeans=False, labels=xdata, positions=xdata) 

    # ...and hacking the SUBBOX plot 
    cap_s = 1.08
    box_s = 1.16
    for i in range(len(ladata[2:])):
        cap = [subxdata[i]/cap_s, subxdata[i]*cap_s]
        subboxplot['caps'][2*i+0].set_data([cap, subboxplot['caps'][2*i+0].get_data()[1]])  #lower cap
        subboxplot['caps'][2*i+1].set_data([cap, subboxplot['caps'][2*i+1].get_data()[1]])  #upper cap

        med = [subxdata[i]/box_s, subxdata[i]*box_s]
        subboxplot['medians'][i].set_data([med, subboxplot['medians'][i].get_data()[1]]) 
        
        box = [subxdata[i]/box_s, subxdata[i]*box_s, subxdata[i]*box_s, subxdata[i]/box_s, subxdata[i]/box_s]
        subboxplot['boxes'][i].set_data([box, subboxplot['boxes'][i].get_data()[1]])

    submedians = [median.get_ydata()[0] for median in subboxplot['medians']]
    ax.scatter(subxdata, submedians, color=colors[draw], marker='o', label='Median', alpha=0.8)  #reduce the marker size! 
    # ...and hacking the BOX plot 
    cap_s = 1.08
    box_s = 1.16
    for i in range(len(ladata)):
        cap = [xdata[i]/cap_s, xdata[i]*cap_s]
        boxplot['caps'][2*i+0].set_data([cap, boxplot['caps'][2*i+0].get_data()[1]])  #lower cap
        boxplot['caps'][2*i+1].set_data([cap, boxplot['caps'][2*i+1].get_data()[1]])  #upper cap

        med = [xdata[i]/box_s, xdata[i]*box_s]
        boxplot['medians'][i].set_data([med, boxplot['medians'][i].get_data()[1]]) 
        
        box = [xdata[i]/box_s, xdata[i]*box_s, xdata[i]*box_s, xdata[i]/box_s, xdata[i]/box_s]
        boxplot['boxes'][i].set_data([box, boxplot['boxes'][i].get_data()[1]])

    medians = [median.get_ydata()[0] for median in boxplot['medians']]
    ax2.scatter(xdata, medians, color=colors[draw], marker='o', label='Median', alpha=0.8)  #reduce the marker size! 

    # for overlaying purposes, such as rhalf noise, rhalf denoised
    if secondfig: 
        medians_two=[]
        submedians_two=[]
        # FOR 15 boxes, all target
        if 1:
            for k in range(len(seconddata)):
                medians_two = np.append(medians_two, np.median(seconddata[k]))
            ax2.scatter(xdata, medians_two, color='r', marker='o', alpha=0.4, label='other median')  #reduce the marker size! 
            for l in range(len(seconddata[2:])):
                submedians_two = np.append(submedians_two, np.median(seconddata[2:][l]))
            ax.scatter(subxdata, submedians_two, color='r', marker='o', alpha=0.4, label='other median')  #reduce the marker size! 
        # FOR 15 boxes, one target - signal only
        if 0:
            medians_two = np.array([np.median(seconddata)]*len(xdata))
            plt.scatter(xdata, medians_two, color='orange', marker='o', alpha=0.5, label='other median')  #reduce the marker size! 

    if thefigure in {'ratios', 'sigmas'}:
        # call the fittings procedure definition
        meds, lows, upps = fitthebox(thefigure, data_lista, theparams)
        # call the fitting definition, which will provide the the_yrecipes below
        ln_ymed, ln_ylow, ln_yupp = plotthefitting(thefigure, meds, lows, upps)
        ymed = np.exp(ln_ymed)
        ylow = np.exp(ln_ylow)
        yupp = np.exp(ln_yupp)
        ax2.plot(the_xrecipe, ymed, color=colors[draw], linestyle='dashed')
        ax2.plot(the_xrecipe, ylow, color=colors[draw], linestyle='dotted', alpha=0.7, label='lower') 
        ax2.plot(the_xrecipe, yupp, color=colors[draw], linestyle='dotted', alpha=0.7, label='upper')  
        ax2.fill_between(the_xrecipe, ylow, yupp, color=colors[draw], alpha=0.2)
    else:
        meds = 0
    
    # SAVING    
    if thefigure == 'rsqrtfull':
        ax.set_yscale('log')
        ax.set_ylabel(r'$R^{1/2}_{noisy}$, $R^{1/2}_{denoised}$')
    if thefigure == 'kvals':
        ax.set_ylabel('k noise')
    if thefigure == 'sigmas':
        ax2.set_yscale('log')
        ax.set_yscale('log')
        ax.set_ylabel(r'$\sigma_{\rho/\rho_o,snr} / \sigma_{\rho/\rho_o,inf}$') 
    if thefigure == 'ratios':
        ax.set_ylabel(r'$R^{1/2}_{denoised}/R^{1/2}_{noisy}$')
    ax.set_xlabel('SNR')

    # really saving now
    plt.savefig('sort_plots/boxplot15_powspec14_%s_i1k_%s'%(thefigure,void_size))
    print('%s_saved!'%thefigure)
    plt.clf()
    return meds
'''