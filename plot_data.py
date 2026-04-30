#!/usr/bin/env python

from globals import *
import matplotlib as mpl
mpl.rcParams.update({
    'errorbar.capsize': 4,
    'lines.markersize': 4
})

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

# ============================= fitting ======================================

def fitthebox(whichfigure, eldata, secondparam=0):
    '''
    GET the median, the 16th ('-' correction), and the 84th ('+' correction) of fit params a, b, c, d, etc
    this then allows us to report param +/- correction 
    '''
    param_a = []
    param_b = []
    param_c = []
    param_d = [] 
    for j in range(15):
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

# ===================================================================

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

# ===================================================================

def errorbarsplot(thefigure, xdata, ydata, void_name=None, secondfig=None, seconddata=None, void_errname=None,
                    styleplot='singular', draw=None, color=None, typeofsigma=None, thefitting=None, void_size=None):
    ax = plt.gca()
    axes = plt.gcf().get_axes()
    if len(axes) > 1:
        ax1 = axes[0]
        ax2 = axes[1]
    if typeofsigma != 'free': ydata = ydata[-1] / ydata # normalise by SNR infinity
    thefigs = [thefigure, secondfig]
    thedata = [ydata, seconddata]
    colores = [color]*4
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
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5, label=None)#labels[num])
                if 1:
                    hfivefile = '../h5files/15cubes15spectra_i1k_1p0_%s_%s.h5'%(void_size, rez) # need to add the powerlaw
                    #hfivefile = '../h5files/15sigma3D_%s.h5'%void_size
                    Fptr = h5py.File(hfivefile,'r')
                    ydata_3D = list(Fptr['sigma_3Ds_%s'%void_size])  
                    ymedians_3D = np.nanpercentile(ydata_3D[0:5], 50)
                    yerrs_3D = [np.nanpercentile(ydata_3D, 16), np.nanpercentile(ydata_3D, 84)]
                    ax.axhline(ymedians_3D, linestyle='dashed', color='b', alpha=0.5, label=r'$\sigma_{\rho/\rho_0}$')
                    xdatafill=np.array([0.006, 1.e-02, 1.e-01, 1.e+00, 3.e+00, 5.e+00, 1.e+01, 1.e+02, 150])
                    ax.fill_between(xdatafill, yerrs_3D[0], yerrs_3D[1], color='blue', alpha=0.1)  #prev alphas: 0.7, 0.2
            else:
                ax.errorbar(xdata[1:], the_medians[1:], yerr=yerrs, linestyle=None, fmt=markers[draw], \
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5, label=labels[num])

        if styleplot =='subpanel':
            if thefigure =='sigmasdenoised' and typeofsigma=='free':
                ax1.errorbar(xdata[1:], the_medians[1:], yerr=yerrs, linestyle=None, fmt=markers[draw], \
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5, label=labels[num])
                if void_name == r'12.5$\mathrm{\mu m}$':
                    ax1.axhline(y=the_medians[k], linestyle='dashed', color=colores[num], alpha=0.5,
                        label=r'$\sigma_{\rho/\rho_0\mathrm{Brunt,denoised}}\mathrm{(SNR\!\to\!\infty)}$')
                else:
                    ax1.axhline(y=the_medians[k], linestyle='dashed', color=colores[num], alpha=0.5)

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
            if 1:
                hfivefile_B = '../h5files/15nonoise_%s.h5'%void_size
                Fptr_B = h5py.File(hfivefile_B, 'r')
                ydata_rhalf = list(Fptr_B['sigma_rhalfs'])
                ymedians_3D = np.nanpercentile(ydata_rhalf, 50)
                yerrs_3D = [np.nanpercentile(ydata_rhalf, 16), np.nanpercentile(ydata_rhalf, 84)]
                xdatafill=np.array([0.06, 1.e-01, 1.e+00, 3.e+00, 5.e+00, 1.e+01, 1.e+02, 150])
            if thefigure =='rsqrtfull':
                print('inside plotting rhalf plot ', void_name)
                ax1.errorbar(xdata[1:], the_medians[1:], yerr=yerrs, linestyle=None, fmt=markers[draw], \
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5, label=labels[num])
                ax1.axhline(ymedians_3D, linestyle='dashed', color=colores[num], alpha=0.5)
                ax1.fill_between(xdatafill, yerrs_3D[0], yerrs_3D[1], color=colores[num], alpha=0.1)
            if thefigure =='rsqrttrue':
                ax2.errorbar(xdata[1:], the_medians[1:], yerr=yerrs, linestyle=None, fmt=markers[draw], \
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5, label=labels[num])
                if void_name == r'12.5$\mathrm{\mu m}$':
                    ax2.axhline(ymedians_3D, linestyle='dashed', color=colores[num], alpha=0.5,
                                label=r'$R^{1/2}\mathrm{(SNR\!\to\!\infty)}$')
                    ax2.fill_between(xdatafill, yerrs_3D[0], yerrs_3D[1], color=colores[num], alpha=0.1)
                else:
                    ax2.axhline(ymedians_3D, linestyle='dashed', color=colores[num], alpha=0.5)
                    ax2.fill_between(xdatafill, yerrs_3D[0], yerrs_3D[1], color=colores[num], alpha=0.1)
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
            meds, lows, upps = fitthebox(thefigure, ydata[1:])
            if thefitting == 'one': 
                medianos1.append(meds)
            if thefitting == 'two': 
                medianos2.append(meds)
            if thefitting == 'three': 
                medianos3.append(meds)

        if void_name == r'100$\mathrm{\mu m}$':
            if thefitting == 'one':
                mediano_1 = np.mean(medianos1,axis=0) 
                ln_ymed, ln_ylow, ln_yupp = plotthefitting(the_xx, thefigure, mediano_1)
            if thefitting == 'two':
                mediano_2 = np.mean(medianos2,axis=0) 
                ln_ymed, ln_ylow, ln_yupp = plotthefitting(the_xx, thefigure, mediano_2)
            if thefitting == 'three':
                if fit =='no':
                    mediano_1 = np.mean(medianos1,axis=0)
                    mediano_2 = np.mean(medianos2,axis=0)
                    mediano = np.append(mediano_1, mediano_2)
                    print('the four params ', mediano)
                    ln_ymed, ln_ylow, ln_yupp = plotthefitting(the_xx, thefigure, mediano)
                if fit =='yes':
                    mediano = np.mean(medianos3,axis=0) 
                    ln_ymed, ln_ylow, ln_yupp = plotthefitting(the_xx, thefigure, mediano)
            ymed = np.exp(ln_ymed)

            if styleplot=='subpanel':
                if whichfigtwo !='sigmasdenoised' and typeofsigma==None:
                    ax2.plot(the_xx, ymed, color='grey', linestyle='dotted', linewidth=3, label=r'$C_1\mathrm{(SNR)}$ fit', alpha=0.7)
                if whichfigtwo=='sigmasdenoised':
                    ax2.plot(the_xx, ymed, color='grey', linestyle='dashed', linewidth=3, label=r'$C_2\mathrm{(SNR)}$ fit', alpha=0.7)
                if typeofsigma=='cthree':
                    if panels=='two':
                        ax2.plot(the_xx, ymed, color='r', linestyle='dashed', linewidth=3, label=r'$C\mathrm{(SNR)}$', alpha=0.7)
                    if panels=='three':
                        ax3.plot(the_xx, ymed, color='r', linestyle='solid', linewidth=3, label=r'$C\mathrm{(SNR)}$', alpha=0.5)
    return  

# ===================================================================

def plotfig(plotstyle, whichfig=None, outfilename=None, typeofsigma=None):
    ax = plt.gca()
    axes = plt.gcf().get_axes()
    if len(axes) > 1:
        ax1 = axes[0]
        ax2 = axes[1]
    if plotstyle == 'singular':
        xlims = (0.006, 150)
        if whichfig == 'kvals':
            label = r'$k_{\mathrm{noise}}/k_{\mathrm{max}}$'
            ylims = (1.5e-2,1.5e0)
            xlims = None
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
        if whichfig == 'sigmas' and typeofsigma == None:
            label = 'Noise correction C(SNR)'
            ax.set_xlabel("SNR", fontsize=12)
            ax.set_ylabel(label, fontsize=12)
            ylims = None
            xlims = None
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
            stop()
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
        ax.set(xscale='log', yscale=yscales, ylim=ylims, xlim=xlims) #EDIT

    if plotstyle == 'subpanel' and whichfig != 'sigmasdenoised':
        ax1.set(xscale='log', yscale='linear') 
        y1label=r'$R^{1/2}$'
        ax1.set_ylabel(y1label, fontsize=21)
        ax1.tick_params(axis='both', which='major', labelsize=21)

        ax2.set(xscale='log', yscale='linear') 
        y2label=r'$R^{1/2}_{\mathrm{denoised}}$'
        ax2.set_ylabel(y2label, fontsize=21)
        ax2.tick_params(axis='both', which='major', labelsize=21)
        ax2.set_xlabel("SNR", fontsize=21)
        ax2.set_xlim(0.06, 150) 

        handles, labels = ax2.get_legend_handles_labels()
        desired_label = r'$R^{1/2}\mathrm{(SNR\!\to\!\infty)}$'
        custom_handle = mlines.Line2D([], [], color='gray', linestyle='dashed', label=desired_label)
        if desired_label in labels:
            handles.append(custom_handle)
            labels.append(desired_label)
            ax2.legend(handles[1:], labels[1:], loc='best', fontsize=20, labelspacing=0.2) 

    if plotstyle == 'subpanel' and whichfig== 'sigmasdenoised':
        ax1.set(xscale='log', yscale='linear') 
        ylims=(0,1)
        ax1.set(ylim=ylims)
        y1label = r'$\sigma_{\mathrm{\rho/\rho_o,Brunt,denoised}}$'
        ax1.set_ylabel(y1label, fontsize=21)
        ax1.tick_params(axis='both', which='major', labelsize=21)

        handles, labels = ax1.get_legend_handles_labels()
        desired_label = r'$\sigma_{\rho/\rho_0\mathrm{Brunt,denoised}}\mathrm{(SNR\!\to\!\infty)}$'
        custom_handle = mlines.Line2D([], [], color='gray', linestyle='dashed', label=desired_label)
        if desired_label in labels:
            #index = labels.index(desired_label)
            #desired_handle = handles.pop(index)
            #desired_label = labels.pop(index)
            handles.append(custom_handle)
            labels.append(desired_label)
            ax1.legend(handles[1:], labels[1:], loc='upper left', fontsize=20, labelspacing=0.2) 

        ax2.set(xscale='log', yscale='log') 
        y2label = r'Noise correction $C_i\mathrm{(SNR)}$'
        ax2.set_ylabel(y2label, fontsize=21)
        ax2.tick_params(axis='both', which='major', labelsize=21)
        ax2.legend(loc='best', fontsize=21)
        if panels=='two':
            ax2.set_xlabel('SNR', fontsize=21)

        if panels=='three':
            ax3.set(xscale='log', yscale='log')
            y3label = r'Noise correction $C\mathrm{(SNR)}$'
            ax3.set_ylabel(y3label, fontsize=21)
            ax3.tick_params(axis='both', which='major', labelsize=21)
            ax3.legend(loc='lower right', fontsize=21)
            ax3.set_xlabel('SNR', fontsize=21)

    # saving Figure
    plt.gcf().tight_layout()
    plt.gcf().savefig(outfilename)
    print(outfilename+' written.', color='magenta')
    plt.clf()


# ===== the following applies in case we are running this in script mode =====
def plot_fig(whichfig, whichfigtwo=None, typeofsigma=None, plotstyle='singular', panels='two', fit = 'no'):

    # figure options:   FIG 3: sigmas (typeofsigma='free', singular)
    #                   FIG 4: kvals2D (typeofsigma=None for all 3 voids, "free" for 50um, 3 rez)
    #                   FIG 5: kvals (singular, typeofsigma=None)
    #                   FIG 6: rsqrtfull (subpanel, two)
    #                   FIG 7: sigmasdenoised (subpanel, three)
    #whichfig = 'kvals2D' # 'kvals', 'kvals2D', 'rsqrtfull', 'rsqrttrue', 'ratios', 'sigmas', 'sigmasdenoised', 'columnddenoise'

    print("=== Plotting "+whichfig+" ===", color='cyan')

    target_snr = [0.01, 0.1, 1, 3, 5, 10, 100, np.inf]
    xdata = np.array(target_snr[:-1])
    thefitting = None

    voids = ['12p5mu', '50mu', '100mu']
    voids_formal = [r'12.5$\mathrm{\mu m}$', r'50$\mathrm{\mu m}$', r'100$\mathrm{\mu m}$']
    voids_eformal = [r'$12.5_{\mu m, \mathrm{denoised}}$', r'$50_{\mu m, \mathrm{denoised}}$', r'$100_{\mu m, \mathrm{denoised}}$']
    markers = ['s', 'o', 'd']
    colors = ['magenta', 'blue', 'green']

    resolutions = ['512', '1024', '2048']
    resolutions = ['512']

    if typeofsigma == 'free':  # or 'free' for the one void figure
        voids_formal = ['512', '1024', '2048']
        voids = ['50mu']
        colors = ['orange', 'blue', 'purple']

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

        for draw, rez in enumerate(resolutions):
            print('resolution(s) = ', rez, color='green')
            for i_void, void_size in enumerate(voids):
                print('void size = ', void_size, color='green')
                hfivefile = '../h5files/15cubes15spectra_i1k_1p0_%s_%s.h5'%(void_size,rez) # need to add the powerlaw
                Fptr = h5py.File(hfivefile, 'r')
                print("read from '"+hfivefile+"'", color='green')

                # define SNR array
                target_snr = np.array([1000000, 0.01, 0.1, 1, 3, 5, 10, 100])

                # get_power spectral data
                kvals = Fptr['kvals2D'][()] # returns a numpy array
                for i in range(len(kvals)):
                    kvals[i] /= kvals[-1] # k/k_max
                P_S = Fptr['thepowspec'][()] # signal spectum
                P_N = Fptr['thenoisepowspec'][()] # noise spectrum
                # getting the k intersections (k_crosses) of S and N
                def get_k_cross(kvals, P_N, P_S):
                    from scipy.interpolate import interp1d
                    k_new = np.linspace(kvals.min(), kvals.max(), 500)
                    f = interp1d(kvals, P_N, kind='linear'); P_N_interp = f(k_new)
                    f = interp1d(kvals, P_S, kind='linear'); P_S_interp = f(k_new)
                    kcross = np.argmin(np.abs(np.log10(P_S_interp/P_N_interp)))
                    P_at_kcross = P_S_interp[kcross]
                    kcross *= kvals[-1] / 500 # scale back to original k grid
                    return kcross, P_at_kcross
                k_crosses_50 = []
                P_at_kcrosses_50 = []
                k_crosses_16 = []
                P_at_kcrosses_16 = []
                k_crosses_84 = []
                P_at_kcrosses_84 = []
                P_N_50 = []
                P_S_50 = []
                P_S_16 = []
                P_S_84 = []
                print("SNR: ", newline=False)
                for i, snr in enumerate(target_snr):
                    print(str(snr)+" ", no_prefix=True, newline=False)
                    P_N_50.append(np.nanpercentile(P_N[i], 50, axis=0))
                    P_S_50.append(np.nanpercentile(P_S[i], 50, axis=0))
                    P_S_16.append(np.nanpercentile(P_S[i], 16, axis=0))
                    P_S_84.append(np.nanpercentile(P_S[i], 84, axis=0))
                    kcross, P_at_kcross = get_k_cross(kvals, P_N_50[-1], P_S_50[0])
                    k_crosses_50.append(kcross)
                    P_at_kcrosses_50.append(P_at_kcross)
                    kcross, P_at_kcross = get_k_cross(kvals, P_N_50[-1], P_S_16[0])
                    k_crosses_16.append(kcross)
                    P_at_kcrosses_16.append(P_at_kcross)
                    kcross, P_at_kcross = get_k_cross(kvals, P_N_50[-1], P_S_84[0])
                    k_crosses_84.append(kcross)
                    P_at_kcrosses_84.append(P_at_kcross)
                print(no_prefix=True)

                # SCATTERS + FITS
                if whichfig != 'kvals2D':

                    # k_noise vs SNR
                    if whichfig == 'kvals':
                        if i_void == 0:
                            fig, ax = plt.subplots(1,1)
                        x = target_snr[2:] # exclude SNR inf and 0.01
                        y = np.array(k_crosses_50[2:])
                        y_lo = y-np.array(np.minimum(k_crosses_16[2:],y)); y_lo[-1] = 0.0
                        y_hi = np.maximum(y,np.array(k_crosses_84[2:]))-y; y_hi[-1] = 0.0
                        ax.errorbar(x, y, yerr=[y_lo, y_hi], linestyle="", marker=markers[i_void], \
                                    color=colors[i_void], alpha=0.6, capsize=2.5, markersize=5, label=voids_formal[i_void])
                        if i_void == 2: # saving condition
                            ax.set(xscale='log', yscale= 'log', xlim=None, ylim=(1.5e-2,1.5e0))
                            ax.set_xlabel(r'SNR', fontsize=12)
                            ax.set_ylabel(r'$k_{\mathrm{noise}}/k_{\mathrm{max}}$', fontsize=12)
                            leg = ax.legend(loc='lower right', fontsize=12, labelspacing=0.2)
                            filename = '../figs/k_noise_%s.pdf'%(rez)
                            plt.savefig(filename)
                            print(filename+' written.', color='magenta')
                            plt.clf()

                    # SINGULAR PLOTS
                    if plotstyle == 'singular':
                        data = np.array(Fptr['%s_%s'%(whichfig,void_size)])
                        print('inside singular')
                        if whichfig == 'rsqrttrue':
                            errorbarsplot(whichfig, xdata, data, void_name=voids_formal[draw], void_errname=voids_eformal[draw], draw=draw, color=colors[draw])
                        if whichfig == 'rsqrtfull':
                            errorbarsplot(whichfig, xdata, data, void_name=voids_formal[draw], void_errname=voids_eformal[draw], draw=draw, color=colors[draw])
                        if whichfig == 'ratios':
                            print('inside fig')
                            theparams = 1  #two params
                            errorbarsplot(whichfig, data, void_name=voids_formal[draw])
                        if whichfig == 'sigmas' and whichfigtwo != 'sigmasdenoised':
                            errorbarsplot(whichfig, xdata, data, void_name=voids_formal[draw], draw=draw, color=colors[draw], typeofsigma='free')
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
                        data = list(Fptr[whichfig])  #returns a list
                        errorbarsplot(whichfig, xdata, data, void_name=voids_formal[draw], void_errname=voids_eformal[draw], styleplot=plotstyle, draw=draw, color=colors[draw], void_size=void_size)
                        print('after rhalf plot ',void_size)

                        whichfig = 'rsqrttrue'
                        data = list(Fptr[whichfig])  #returns a list
                        errorbarsplot(whichfig, xdata, data, void_name=voids_formal[draw], void_errname=voids_eformal[draw], styleplot=plotstyle, draw=draw, color=colors[draw], void_size=void_size)
                        #              secondfig=whichfigtwo, seconddata=figtwodata)
                        print('after true plot ', void_size)
                        whichfig = 'rsqrtfull'

                if whichfig not in ['kvals', 'kvals2D']:
                    outfilename = '../figs/'+whichfig+'_'+void_size+'_'+str(rez)+'.pdf'
                    plotfig(plotstyle, whichfig=whichfig, outfilename=outfilename, typeofsigma=typeofsigma)

                # POWERSPECTRA PLOTS
                if whichfig == 'kvals2D':
                    # prep fig and read spectral data
                    fig, ax = plt.subplots(1,1)
                    # plot S spectrum
                    ax.plot(kvals[1:], P_S_50[0][1:], color='r', linewidth=0.8, label='S', zorder=9)
                    ax.fill_between(kvals[1:], P_S_16[0][1:], P_S_84[0][1:], color='red', alpha=0.2, edgecolor='none', zorder=9)
                    # loop over SNRs
                    for i, snr in enumerate(target_snr):
                        if snr <= 100:
                            label = None
                            # plot N spectra and SNR labels
                            if snr == 0.01: label = 'N'
                            ax.plot(kvals[1:], P_N_50[i][1:], color='grey', linewidth=0.8, linestyle='dashed', label=label)
                            yoff = 0.0
                            if snr == 5: yoff = -5e-4
                            plt.text(kvals[-1]+0.05, P_N_50[i][-1]+yoff, str(snr), fontsize=11, ha='left', va='center')
                            # plot S+N spectra
                            if snr == 0.01: label = 'S+N'
                            ax.plot(kvals[1:], P_S_50[i][1:], color=colors[i_void], linewidth=0.8, alpha=0.7, label=label)
                        if snr <= 0.01 or snr >= 100: continue # skip SNRs that don't have useful k crossings
                        ax.plot([k_crosses_50[i],k_crosses_50[i]], [0,P_at_kcrosses_50[i]], color=colors[i_void],
                                    linestyle=(0, (1, 1)), alpha=0.8, zorder=10)
                        ax.fill_between([k_crosses_16[i], max([k_crosses_50[i], k_crosses_84[i]])], [0, 0],
                                            [P_at_kcrosses_50[i], P_at_kcrosses_50[i]],
                                            color=colors[i_void], alpha=0.1, edgecolor='none', zorder=10)
                    # final plot settings and saving
                    ax.text(0.905, 0.94, 'SNR', transform=ax.transAxes)
                    ax.set(xscale='log', yscale= 'log', xlim=(2e-2,1.6), ylim=(1e-7,2e4))
                    ax.set_xlabel(r'$k/k_{\mathrm{max}}$', fontsize=12)
                    ax.set_ylabel(r'$P_{\mathrm{2D}}$', fontsize=12)
                    if void_size == '100mu':
                        bbox_to_anchor=(0.01, 1.02)
                    else:
                        bbox_to_anchor=(0.01, 1.01)
                    main_legend = ax.legend(loc='upper left', ncol=1, fontsize=12, labelspacing=0.1, bbox_to_anchor=bbox_to_anchor)
                    ax.add_artist(main_legend)
                    custom_line = mlines.Line2D([], [], color=colors[i_void], linestyle=(0, (1, 1)), alpha=0.8, label=r'$k_{\mathrm{noise}}/k_{\mathrm{max}}$')
                    k_legend = ax.legend(handles=[custom_line], loc='upper left', fontsize=12, bbox_to_anchor=(0.23,1.01))
                    filename = '../figs/spectra_%s_%s.pdf'%(void_size, rez)
                    plt.savefig(filename)
                    print(filename+' written.', color='magenta')
                    plt.clf()


    # TABLE 2, and can find TABLE 1 values.
    if 0:
        voids = ['12p5mu', '50mu', '100mu']
        for num, void_size in enumerate(voids):
            hfivefile = '../h5files/15sigma3D_%s.h5'%void_size
            hfivefile_B = '../h5files/15nonoise_%s.h5'%void_size
            Fptr = h5py.File(hfivefile,'r')
            Fptr_B = h5py.File(hfivefile_B,'r')
            ydata_3D = list(Fptr['sigma_3Ds_%s'%void_size])
            ydata_Brunt = list(Fptr_B['sigma_brunts_%s'%void_size])
            ydata_2Ds = list(Fptr_B['sigma_2Ds_%s'%void_size])
            ydata_rhalf = list(Fptr_B['sigma_rhalfs_%s'%void_size])
            fourlists = [ydata_3D, ydata_2Ds, ydata_rhalf, ydata_Brunt]
            median_vals =[]
            lower_vals =[]
            upper_vals =[]
            for num, param in enumerate(fourlists):
                median_vals.append(np.percentile(fourlists[num],50))
                lower_vals.append(np.percentile(fourlists[num],16))
                upper_vals.append(np.percentile(fourlists[num],84))
            print('void size is ',void_size)
            print('median values ',np.array(median_vals))
            print('upper val error ',np.array(upper_vals)-np.array(median_vals))
            print('lower val error ',np.array(median_vals)-np.array(lower_vals))

            # Section 3.2
            med1 = np.percentile([ydata_Brunt[0],ydata_Brunt[5],ydata_Brunt[10]],50)
            med2 = np.percentile([ydata_Brunt[1],ydata_Brunt[6],ydata_Brunt[11]],50)
            med3 = np.percentile([ydata_Brunt[2],ydata_Brunt[7],ydata_Brunt[12]],50)
            med4 = np.percentile([ydata_Brunt[3],ydata_Brunt[8],ydata_Brunt[13]],50)
            med5 = np.percentile([ydata_Brunt[4],ydata_Brunt[9],ydata_Brunt[14]],50)
            print('medians ',med1, med2, med3, med4, med5)
            medstd = np.std([med1,med2, med3, med4, med5])

            # Section 3.3 - relative errors sigr vs sigr_b
            rel_error = abs(np.array(ydata_Brunt)-np.array(ydata_3D))/abs(np.array(ydata_3D))
            rel_error_percent = rel_error*100 
            rel_error_mean = np.mean(rel_error_percent)
            rel_error_std = np.std(rel_error_percent)

    # APPENDIX FIGURES
    if 0:
        void_size = '50mu'
        hfivefile_0p0 = '../h5files/15cubes15spectra_i1k_0p0_%s.h5'%void_size
        hfivefile_0p5 = '../h5files/15cubes15spectra_i1k_0p5_%s.h5'%void_size
        hfivefile_1p0 = '../h5files/15cubes15spectra_i1k_1p0_%s.h5'%void_size
        hfivefile_1p5 = '../h5files/15cubes15spectra_i1k_1p5_%s.h5'%void_size
        hfivefile_2p0 = '../h5files/15cubes15spectra_i1k_2p0_%s.h5'%void_size
        hfivefiles = [hfivefile_0p0, hfivefile_0p5, hfivefile_1p0, hfivefile_1p5, hfivefile_2p0]
        the_medians = []
        the_uppers = []
        the_lowers = []
        target_pl = [0.0, 0.5, 1.0, 1.5, 2.0]
        pl_data = np.array(target_pl)
        target_snr = [0.1, 1, 3, 5, 10, 100]
        colors=['k','gray','b','g','orange','r']
        collect_snr = defaultdict(list)
        collect_lowersnr= defaultdict(list)
        collect_uppersnr = defaultdict(list)

    # TABLES VALUES EXPLORATION
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
                hfivefile = '../h5files/15sigma3D_%s.h5'%void_size
                hfivefile_B = '../h5files/15nonoise_%s.h5'%void_size
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
            fig.savefig('../figs/15boxes_brunt_3D_%s'%void_size)
            plt.close()
            Fptr.close()

    # FIGURE 7 EXPLORATION with noise spec k=1, k=1/2, k=0, k=2, one void size only
    # WARNING: changed lines 162 - 166
    if 0:
        voids_formal = ['pl=1', 'pl=1/2', 'pl=0', 'pl=2']
        plotstyle='subpanel'
        whichfig = 'sigmasdenoised' 
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



def plot_noisy_map_examples(turb_data):
    # plot example maps of noisy data
    for isnr in [2,3,4,6]:
        # check SNR match
        dat = turb_data['data'][isnr][1] # extract target SNR data for centre box and y LOS (index 1)
        vmin, vmax = 5e-3, 1e-2
        cmap = 'magma'
        SNR = SNR_list[isnr]
        if dat.SNR[0] != SNR: stop()
        save = figpath+'noisy_data_map_SNR_'+cfp.eform(SNR,prec=1)+'_'+void_size+'.pdf'
        po = cfp.plot_map(dat.noisy_data[0], axes_format=["", ""], axes_ticks=[[],[]],
                            vmin=vmin, vmax=vmax, cmap=cmap, colorbar=False)
        po.plt().text(0.5, 0.95, r"$\mathrm{SNR} = "+str(SNR)+r"$", ha='center', va='center', transform=po.plt().gca().transAxes,
                        bbox=dict(boxstyle="round", fc="white", alpha=0.8))
        cfp.show_or_save_plot(save=save)
    # plot colorbar
    cfp.plot_colorbar(vmin=vmin, vmax=vmax, cmap=cmap, label=r'Column density $(\mathrm{g\,cm^{-2}})$',
                        panels=3, orientation="horizontal", save=figpath+'noisy_data_map_cbar_'+void_size+'.pdf')


def plot_sigrho_Brunt_vs_SNR(turb_data):
    # figure 3 (Brunt sigma_rho as function of SNR without noise correction)
    for isnr, SNR in enumerate(SNR_list):
        x = SNR
        y = turb_data['stats'][isnr].avg.Brunt_sigma_3D
        yerr = [y-turb_data['stats'][isnr].p16.Brunt_sigma_3D, turb_data['stats'][isnr].p84.Brunt_sigma_3D-y]
        cfp.plot(x=x, y=y, yerr=np.array([yerr]).T, color=color, marker='o', alpha=0.8)
    # draw true sigma_rho
    xlim = [5e-3,2e2]
    y = turb_data['stats'][0].avg.sigma_3D
    yerr = [y-turb_data['stats'][0].p16.sigma_3D, turb_data['stats'][0].p84.sigma_3D-y]
    cfp.plot(x=xlim, y=[y,y], yerr=np.array([yerr,yerr]).T, shaded_err=True,
                linestyle='--', color=color, label=r'True $\sigma_{\rho/\rho_0}$')
    # generate final plot and save
    cfp.plot(xlabel=r"$\mathrm{SNR}$", ylabel=r"$\sigma_{\rho/\rho_0,\mathrm{Brunt}}$",
                xlog=True, ylog=True, xlim=xlim, ylim=[1e-1,2e2], legend_loc='upper right',
                save=figpath+'Brunt_sigrho_vs_SNR_'+void_size+'.pdf')


# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    res = 1024 # resolution
    void_sizes = ['12p5mu', '50mu', '100mu']
    void_sizes = ['50mu']
    void_size = void_sizes[0]
    color = 'tab:blue'

    # load data
    dillfile = get_filename(res, void_size, dill=True)
    with open(dillfile, "rb") as f:
        turb_data = dill.load(f)
    print("Read '"+dillfile+"'.", color='green')

    # figure 3 (Brunt sigma_rho as function of SNR without noise correction)
    plot_sigrho_Brunt_vs_SNR(turb_data)

    # figure 2 (noisy map examples)
    plot_noisy_map_examples(turb_data)

    stop()

    # figure options:   FIG 3: sigmas (typeofsigma='free', singular)
    #                   FIG 4: kvals2D (typeofsigma=None for all 3 voids, "free" for 50um, 3 rez)
    #                   FIG 5: kvals (singular, typeofsigma=None)
    #                   FIG 6: rsqrtfull (subpanel, two)
    #                   FIG 7: sigmasdenoised (subpanel, three)

    # args: 'kvals', 'kvals2D', 'rsqrtfull', 'rsqrttrue', 'ratios', 'sigmas', 'sigmasdenoised', 'columnddenoise'

    plot_fig('rsqrtfull', plotstyle='subpanel', panels='two') # Fig. 6
    plot_fig('kvals2D')
    plot_fig('kvals') # Fig. 5
