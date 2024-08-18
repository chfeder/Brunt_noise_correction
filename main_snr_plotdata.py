
from startlaser import * 
from collections import defaultdict
import matplotlib.ticker as ticker
#from labelLines import *
# GET INFO FROM STORAGE 


# POSSIBLE FIT FUNCTIONS
def linfunct(x, a, b):  
    y = a * x + b 
    #y = np.log(pre_y)
    return y

def lininvfunct(x, a, b=1):
    #pre_y = a * x**(-1) + b
    y = a * x**(-1) + b
    #y = np.log(pre_y)
    return y

def tanfunct_simple(x, first, second):  #used to be x, c, d, e; try this first when data is ready
    #y = first * np.tanh(second*np.log10(x)) + (1-first)  #c should be negative, d positive; c+e=1
    y = first * np.tanh(second*x) + (1-first)  #c should be negative, d positive; c+e=1
    return y

def gaussianlike(x, a, b): 
    #y = a * np.exp(-(x-1)**2 / (2 * b**2)) + c
    y = a * np.exp(-(x)**2 /b) + 1
    return y

def sigmarho_fit_old_old(x,a):  
    #y = 1/ (lininvfunct(x, a, b=1.0) * tanfunct_simple(x, first=first, second=second)) 
    #y = (lininvfunct(x, a, b=1) * tanfunct_simple(x, first=first, second=second)) 
    pre_y = (lininvfunct(x, a, b=1) * tanfunct_simple(x, first=first, second=second)) 
    y = np.log(pre_y)
    return y 

def sigmarho_fit_old(x, a, b, c):
    y = 1/(a * x**(b)+c)
    return y 

def sigmarho_fit(x, a, b):
    y = np.log( x**a / (x**a + b) )
    #y = (x**a) / (x**a + b)
    return y 



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
    param_e = []
    for j in range(numofboxes):  
        if whichfigure == 'sigmas':
            xdata = [0.1, 1, 3, 5, 10, 100] 
        oneboxproj = []
        for k in range(len(xdata)):
            prexdata = eldata[k][j]
            oneboxproj = np.append(oneboxproj, prexdata) 
        
        if whichfigure == 'ratios': 
            #mask = ~np.isnan(oneboxproj)
            #xdata_filtered = np.array(xdata)[mask]
            #ydata_filtered = oneboxproj[mask]
            print('inside the ratios')
            secondparam = "ratios"
            #params = {"first":[-2,-1,-0.1], "second":[0.5,2.0,4]} 
            rets = cfp.fit(tanfunct_simple, xdata[2:], oneboxproj[2:])#, params=params) 
            #rets = cfp.fit(tanfunct_simple, xdata[2:], oneboxproj[2:])#, params=params) 

        if whichfigure == 'sigmas':
            print('inside the sigmas')
            secondparam = "ratios"
            #params= {'a':[2, 3, 4], 'b':[-2,-1,0], 'c':[0.5, 0.75, 1], "first":[-2,-1,-0.1], "second":[0.5,2.0,4]}
            #params= {'a':[1.25, 1.5, 1.75], 'b':[-1.35,-1.25,-1.15]}#, 'c':[0.15,0.0,0.25], 'd':[0.05,0.0,0.15]}
            #params= {'a':[1.4, 1.6, 1.8], 'b':[-1.35,-1.15,-0.95]}
            #params=None
            params = {"a": [0.1, 1, 3.0], "b": [0.1, 1, 3.0]}
            #params= {'a':[-1, -0.5, 0], 'b':[0.5,0.8,1], 'c':[-1, -0.5, -0.01]}
            rets = cfp.fit(sigmarho_fit, xdata, np.log(oneboxproj), params=params)

        if whichfigure == 'kvals':
            print('inside fitting for kvals')
            secondparam = 1
            rets = cfp.fit(linfunct, xdata[1:-1], oneboxproj[1:-1])
            print('how many params ',secondparam)

        param_a = np.append(param_a, rets.popt[0])
        if secondparam == 'gauss' or secondparam == 'ratios':
            param_b = np.append(param_b, rets.popt[1])
        if secondparam== 'gauss':
            param_c = np.append(param_c, rets.popt[2])
            #param_d = np.append(param_d, rets.popt[3])
            #param_e = np.append(param_e, rets.popt[4])
            #param_f = np.append(param_f, rets.popt[5])

    med_a = np.percentile(param_a, 50) #same as np.median()
    p16_a = np.percentile(param_a, 16)
    p84_a = np.percentile(param_a, 84)
    if secondparam == 'gauss' or secondparam =='ratios':
        med_b = np.percentile(param_b, 50)
        p16_b = np.percentile(param_b, 16)
        p84_b = np.percentile(param_b, 84)
    if secondparam == 'gauss':
        med_c = np.percentile(param_c, 50)
        p16_c = np.percentile(param_c, 16)
        p84_c = np.percentile(param_c, 84)
        #med_d = np.percentile(param_d, 50)
        #p16_d = np.percentile(param_d, 16)
        #p84_d = np.percentile(param_d, 84)
        #med_e = np.percentile(param_e, 50)
        #p16_e = np.percentile(param_e, 16)
        #p84_e = np.percentile(param_e, 84)
        #med_f = np.percentile(param_f, 50)
        #p16_f = np.percentile(param_f, 16)
        #p84_f = np.percentile(param_f, 84)

    if secondparam == 0:
        med_b = p16_b = p84_b = 0
    if 0: 
        med_c = p16_c = p84_c = 0
    if 1:
        med_d = p16_d = p84_d = 0
    if 1:
        med_e = p16_e = p84_e = 0

    if secondparam == 0:
        param_meds = [med_a]
        param_16s = [p16_a] 
        param_84s = [p84_a] 
    else:
        print("there is a b and or c")
        param_meds = [med_a, med_b]#, med_c]#, med_d]#, med_e]#, med_f]
        param_16s = [p16_a, p16_b]#, p16_c]#, p16_d]#, p16_e]#, p16_f]  #report these as - for param
        param_84s = [p84_a, p84_b]#, p84_c]#, p84_d]#, p84_e]#, p84_f]  #report these as + for param
    
    return param_meds, param_16s, param_84s


# TO PLOT THE FITTING 
def plotthefitting(the_xxx, thefigures, themedians, sixteenths, eightyfours):
    if thefigures == 'sigmas':
        the_yrecipemed = sigmarho_fit(the_xxx, *themedians)  
        the_yrecipelow = sigmarho_fit(the_xxx, *sixteenths)
        the_yrecipeupp = sigmarho_fit(the_xxx, *eightyfours)
    if thefigures == 'ratios':
        the_yrecipemed = tanfunct_simple(the_xrecipe, *themedians)  
        the_yrecipelow = tanfunct_simple(the_xrecipe, *sixteenths)
        the_yrecipeupp = tanfunct_simple(the_xrecipe, *eightyfours)
    if thefigures == 'kvals':
        the_yrecipemed = linfunct(the_xrecipe, *themedians)  
        the_yrecipelow = linfunct(the_xrecipe, *sixteenths)
        the_yrecipeupp = linfunct(the_xrecipe, *eightyfours)
    return the_yrecipemed, the_yrecipelow, the_yrecipeupp


# MAKE A BOXPLOT AND SCATTER PLOT 
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
    subboxplot = ax.boxplot(data_lista[2:], showfliers=True, showmeans=False, labels=subxdata, positions=subxdata) #\
                        #, boxprops=dict(color="none"), capprops=dict(color="none"), medianprops=dict(color="none"))  
    boxplot = ax2.boxplot(data_lista, showfliers=True, showmeans=False, labels=xdata, positions=xdata) #\
                        #, boxprops=dict(color="none"), capprops=dict(color="none"), medianprops=dict(color="none"))  


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

        # to be EDITED; same as above
        #stop()
        '''
        meds, lows, upps = fitthebox(thefigure, data_lista, theparams)
        ymed, ylow, yupp = plotthefitting(thefigure, meds, lows, upps)
        ax.plot(the_xrecipe, ymed, color=colors[draw], linestyle='dashed')
        ax.plot(the_xrecipe, ylow, color=colors[draw], linestyle='dotted', alpha=0.7, label='lower') 
        ax.plot(the_xrecipe, yupp, color=colors[draw], linestyle='dotted', alpha=0.7, label='upper')  
        ax.fill_between(the_xrecipe, ylow, yupp, color=colors[draw], alpha=0.2)
        '''
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


def boxplotscatterone(thefigure, ladata, theparams=0, secondfig=None, seconddata=None):
    fig,ax = plt.subplots(1,1)
    ax.set_xscale('log')
    ax.set_xlim(np.array(xdata).min()/1.20, np.array(xdata).max()*1.20)

    if thefigure == 'sigmas':  
        print("inside sigmas!!")
        the_data = ladata/ladata[-1]
        data_lista = the_data.tolist() 
    elif thefigure == 'ratios':
        print("inside ratios")
        pre_data_lista = ladata[1:]
    else:
        print("inside else")
        data_lista = ladata
    boxplot = ax.boxplot(data_lista, showfliers=True, showmeans=False, labels=xdata, positions=xdata)

    cap_s = 1.08
    box_s = 1.16
    for i in range(len(data_lista)):
        cap = [xdata[i]/cap_s, xdata[i]*cap_s]
        boxplot['caps'][2*i+0].set_data([cap, boxplot['caps'][2*i+0].get_data()[1]])  #lower cap
        boxplot['caps'][2*i+1].set_data([cap, boxplot['caps'][2*i+1].get_data()[1]])  #upper cap

        med = [xdata[i]/box_s, xdata[i]*box_s]
        boxplot['medians'][i].set_data([med, boxplot['medians'][i].get_data()[1]]) 
        
        box = [xdata[i]/box_s, xdata[i]*box_s, xdata[i]*box_s, xdata[i]/box_s, xdata[i]/box_s]
        boxplot['boxes'][i].set_data([box, boxplot['boxes'][i].get_data()[1]])

    medians = [median.get_ydata()[0] for median in boxplot['medians']]
    ax.scatter(xdata, medians, color=colors[draw], marker='o', label='Median', alpha=0.8)  #reduce the marker size! 


    # for overlaying purposes, such as rhalf noise, rhalf denoised
    if secondfig: 
        medians_two=[]
        submedians_two=[]
        # FOR 15 boxes, all target
        if 1:
            for k in range(len(seconddata)):
                medians_two = np.append(medians_two, np.nanmedian(seconddata[k]))
            ax.scatter(xdata, medians_two, color='r', marker='o', alpha=0.4, label='other median')  #reduce the marker size! 
            #stop()

    if thefigure in {'sigmas'}:#, 'kvals'}:
        # call the fittings procedure definition
        meds, lows, upps = fitthebox(thefigure, data_lista, theparams)
        # call the fitting definition, which will provide the the_yrecipes below
        if thefigure == 'sigmas':
            ln_ymed, ln_ylow, ln_yupp = plotthefitting(thefigure, meds, lows, upps)
            ymed = np.exp(ln_ymed)
            ylow = np.exp(ln_ylow)
            yupp = np.exp(ln_yupp)
            the_xrecipe = np.linspace(xdata.min(), xdata.max(), num=2500)  
        elif thefigure == 'kvals':
            ymed, ylow, yupp = plotthefitting(thefigure, meds, lows, upps)
            the_xrecipe = np.linspace(xdata[1:-1].min(), xdata[1:-1].max(), num=2500)  
        #stop()
        ax.plot(the_xrecipe, ymed, color=colors[draw], linestyle='dashed')
        ax.plot(the_xrecipe, ylow, color=colors[draw], linestyle='dotted', alpha=0.7, label='lower') 
        ax.plot(the_xrecipe, yupp, color=colors[draw], linestyle='dotted', alpha=0.7, label='upper')  
        ax.fill_between(the_xrecipe, ylow, yupp, color=colors[draw], alpha=0.2)
    else:
        meds = 0
    
    # SAVING    
    if thefigure == 'rsqrttrue':
        ax.set_yscale('log')
        ax.set_ylabel(r'$R^{1/2}_{denoised}$')
    if thefigure == 'rsqrtfull':
        ax.set_yscale('log')
        ax.set_ylabel(r'$R^{1/2}_{noisy}$, $R^{1/2}_{denoised}$')
        #ax.set_ylabel(r'$R^{1/2}_{noisy}$')
    if thefigure == 'kvals':
        #ax.set_ylabel('k noise')
        ax.set_ylabel(r'$k_{noise}/k_{max}$')
        ax.set_yscale('log')
        #ax.set_ylim(1e-1,1e2)
    if thefigure == 'sigmas':
        ax.set_yscale('log')
        ax.set_ylabel(r'$\sigma_{\rho/\rho_o,snr} / \sigma_{\rho/\rho_o,inf}$') 
        #ax.set_ylabel(r'$\sigma_{\/rho/\rho_o,snr} / \sigma_{\rho/\rho_o,inf}$') 
    if thefigure == 'ratios':
        ax.set_ylabel(r'$R^{1/2}_{denoised}/R^{1/2}_{noisy}$')
        ax.set_yscale('log')
    ax.set_xlabel('SNR')

    # really saving now
    if 1:  #quick brunt test with real data
        plt.savefig('sort_plots/boxplot15_powspec14_%s_i1k_%s'%(thefigure,void_size))
        print('%s_saved!'%thefigure)
        plt.clf()
    return meds


def errorbarsplot(thefigure, ladata, void_name=None, secondfig=None, seconddata=None, void_errname=None, styleplot='singular'):
    thefigs = [thefigure, secondfig]
    thedata = [ladata, seconddata]
    colores = [colors[draw], colors[draw]]
    labels = [void_name, void_errname]
    alphas = [0.7, 0.3]
    markers = ['s', 'o', 'D']
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
            if typeofsigma == 'free':
                ax.errorbar(xdata, the_medians, yerr=[lower_error, upper_error], linestyle=None, fmt=markers[1], \
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5, label=labels[num])
            else:
                ax.errorbar(xdata[1:], the_medians[1:], yerr=yerrs, linestyle=None, fmt=markers[draw], \
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5, label=labels[num])
        if styleplot =='subpanel':
            if thefigure =='rsqrtfull':
                ax1.errorbar(xdata[1:], the_medians[1:], yerr=yerrs, linestyle=None, fmt=markers[draw], \
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5, label=labels[num])
                #ax1.axhline(y=the_medians[k], linestyle='dashed', color=colores[num], alpha=0.4)
            if thefigure =='rsqrttrue':
                ax2.errorbar(xdata[1:], the_medians[1:], yerr=yerrs, linestyle=None, fmt=markers[draw], \
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5, label=labels[num])
                #ax2.axhline(y=the_medians[k], linestyle='dashed', color=colores[num], alpha=0.4)
            if thefigure =='ratios':
                ax3.errorbar(xdata[1:], the_medians[1:], yerr=yerrs, linestyle=None, fmt=markers[draw], \
                            ecolor=colores[num], color=colores[num], alpha=alphas[num], capsize=5, label=labels[num])
                if void_name =='12.5$\mathrm{\mu m}$':
                    ax3.axhline(y=1, linestyle='dashed', color='gray', alpha=0.7, label=r'$R^{1/2}_{\mathrm{denoised}}\mathrm{(SNR\!\to\!\infty)}$')
        if secondfig == None:
            break
        if secondfig:
            break

    if thefigure == 'sigmas' and typeofsigma==None:
        if 0:  #fitting just the medians of 15 boxes
            params = {"a": [0.1, 1, 3.0], "b": [0.1, 1, 3.0]}
            #params = None
            rets = cfp.fit(sigmarho_fit, xdata, np.log(the_medians), params=params)
            xfit = np.linspace(xdata.min(), xdata.max(), num=500)
            yfit = np.exp(sigmarho_fit(xfit, *rets.popt))
            ax.plot(xfit, yfit, color=colors[draw], linestyle='dashed', label=labels[num], alpha=0.5)

        if 1:  #detailed, fitting each box, then...
            the_xx = np.linspace(xdata[1:].min(), xdata[1:].max(), num=2500)  
            # call the fittings procedure definition
            meds, lows, upps = fitthebox(thefigure, ladata[1:])
            print('median a, b, and maybe c', meds)
            # call the fitting definition, which will provide the the_yrecipes below

            mediano = np.append(medianos, meds) 
            arribas = np.append(arriba, upps)
            abajos = np.append(abajo,lows)
            if void_name == '100$\mathrm{\mu m}$':
                ln_ymed, ln_ylow, ln_yupp = plotthefitting(the_xx, thefigure, mediano, abajos, arribas)
                #ln_ymed, ln_ylow, ln_yupp = plotthefitting(the_xx, thefigure, meds, lows, upps)
                ymed = np.exp(ln_ymed)
                #ylow = np.exp(ln_ylow)
                #yupp = np.exp(ln_yupp)
                ax.plot(the_xx, ymed, color='k', linestyle='dashed', linewidth=3, label='fit', alpha=0.5)
                #ax.plot(the_xx, ymed, color=colors[draw], linestyle='dashed', label=labels[num], alpha=0.5)

                #ax.plot(the_xrecipe, ln_ylow, color=colors[draw], linestyle='dotted', alpha=0.7, label='lower') 
                #ax.plot(the_xrecipe, ln_yupp, color=colors[draw], linestyle='dotted', alpha=0.7, label='upper')  
                #ax.fill_between(the_xrecipe, ln_ylow, ln_yupp, color=colors[draw], alpha=0.2)

    return  



def plotfig(plotstyle):
    if plotstyle == 'singular':
        if whichfig == 'kvals':
            label = r'$k_{\mathrm{noise}}/k_{\mathrm{max}}$'
            ylims = None#(1.5e-2,1.5e0)
            yscales = 'log'
            ax.set_xlabel("SNR", fontsize=13)
            ax.set_ylabel(label, fontsize=13)
            ax.legend(loc='lower right', fontsize=12, labelspacing=0.2)
        if whichfig == 'rsqrtfull':
            label = r'$R^{1/2}_{\mathrm{noisy}}$'
            ylims = None
            yscales = 'linear'
        if whichfig == 'rsqrttrue':
            label = r'$R^{1/2}_{\mathrm{denoised}}$'
            ylims = None
            yscales = 'linear'
        if whichfig == 'ratios':
            label = r'$R^{1/2}_{\mathrm{denoised}}/R^{1/2}_{\mathrm{noisy}}$'
            ylims = None
            yscales = 'linear'#'log'
        if whichfig == 'sigmas':
            #label = r'$C\mathrm{(SNR)}=\sigma_{\mathrm{\rho/\rho_o,3D}} (\mathrm{SNR}\!\to\!\infty)  / \sigma_{\mathrm{\rho/\rho_o,3D}}\mathrm{(SNR)}$'  
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
        if typeofsigma == 'free':
            label = r'$\sigma_{\mathrm{\rho/\rho_o,\,3D}}$' #\, is a thin space, \: a medium space  
            ax.set_xlabel("SNR", fontsize=13)
            ax.set_ylabel(label, fontsize=13)
            ax.legend(loc='best', fontsize=12)

        ax.set(xscale='log', yscale=yscales, ylim=ylims)
        #ax.set(xscale='log', yscale=yscales, ylabel=label, xlabel='SNR', ylim=ylims)

        #if whichfig != 'sigmas':
        #    ax.legend(loc='best', fontsize=14, labelspacing=0.2)

    if plotstyle == 'subpanel':
        ax1.set(xscale='log', yscale='linear') 
        y1label=r'$R^{1/2}_{\mathrm{noisy}}$'
        ax1.set_ylabel(y1label, fontsize=21)

        ax2.set(xscale='log', yscale='linear') 
        y2label=r'$R^{1/2}_{\mathrm{denoised}}$'
        ax2.set_ylabel(y2label, fontsize=21)

        ax3.set(xscale='log', yscale='linear')
        y3label=r'$R^{1/2}_{\mathrm{denoised}}/R^{1/2}_{\mathrm{noisy}}$'
        ax3.set_ylabel(y3label, fontsize=21)
        ax3.set_xlabel("SNR", fontsize=21)

        handles, labels = ax3.get_legend_handles_labels()
        desired_label = r'$R^{1/2}_{\mathrm{denoised}}\mathrm{(SNR\!\to\!\infty)}$'
        if desired_label in labels:
            index = labels.index(desired_label)
            desired_handle = handles.pop(index)
            desired_label = labels.pop(index)
            handles.append(desired_handle)
            labels.append(desired_label)
            ax3.legend(handles, labels, loc='best', fontsize=20, labelspacing=0.2) 


    # really saving now
    fig.tight_layout()   #trouble for sigmas?
    fig.savefig('sort_plots/boxplot15_powspec14_%s_i1knobox.pdf'%(whichfig))
    print('%s_saved!'%whichfig)
    plt.clf()




# MAIN: LOOP FOR THREE VOIDS!
numofboxes = 15
#subxdata = np.array(target_snr[2:-1])
#the_subxrecipe = np.linspace(subxdata.min(), subxdata.max(), num=2500)  

whichfig = 'kvals2D'
typeofsigma = None#'free'
target_snr = [0.01, 0.1, 1, 3, 5, 10, 100, inf]  
#if typeofsigma == 'free' or whichfig=='kvals2D':
#    target_snr = [0.01, 0.1, 1, 3, 5, 10, 100, inf]  
#else:
#    target_snr = [0.1, 1, 3, 5, 10, 100, inf]  
xdata = np.array(target_snr[:-1])
the_xrecipe = np.linspace(xdata.min(), xdata.max(), num=2500)  

k_plot = 'divkmax'  # divkmax or knoise
whichfigtwo = None#'rsqrttrue'
mainfigure = None#'sigmas'
plotstyle = 'singular'  # singular or subpanel
boxorplot = 'plot'  #dismissing the boxplots...

voids = ['12p5mu', '50mu', '100mu']
voids_formal = ['12.5$\mathrm{\mu m}$', '50$\mathrm{\mu m}$', '100$\mathrm{\mu m}$']
voids_eformal = ['$12.5_{\mu m, \mathrm{denoised}}$', '$50_{\mu m, \mathrm{denoised}}$', '$100_{\mu m, \mathrm{denoised}}$',]
colors = ['magenta', 'b', 'g']

if typeofsigma == 'free':  #or for debug 
    voids = ['50mu']  
    voids_formal = ['50$\mathrm{\mu m}$']
    colors = ['b']    

if 1:
    if whichfig != 'kvals2D' and plotstyle=='singular':
        fig,ax = plt.subplots(1,1)  #TO OVERLAY 
    if whichfig != 'kvals2D' and plotstyle=='subpanel':
        fig,(ax1, ax2, ax3) = plt.subplots(3,1, sharex=True, figsize=(8, 12))  
    if whichfig == 'sigmas' and typeofsigma == None:
        medianos = []
        arriba = []
        abajo = []
    for draw, void_size in enumerate(voids):
        print('void size is ',void_size)
        hfivefile = '/home/552/lv9897/sflaser/h5files/boxplot15_powspec14_i1k_SintNstartk0_%s.h5'%void_size
        Fptr = h5py.File(hfivefile,'r')

        # ALL OTHER DATA: BOXPLOTS + SCATTERS + FITS
        if whichfig != 'kvals2D':

            # SINGULAR
            if plotstyle == 'singular':
                data = list(Fptr['%s_%s'%(whichfig,void_size)])  #returns a list 
                print('inside singular')

                if whichfig == 'kvals' and k_plot=='divkmax':
                    for val, kvalues in enumerate(data):
                        data[val] = data[val]/data[-1] 
                    errorbarsplot(whichfig, data, void_name=voids_formal[draw])  
                if whichfig == 'rsqrttrue':
                    if boxorplot == 'plot':
                        errorbarsplot(whichfig, data, void_name=voids_formal[draw], void_errname=voids_eformal[draw])  
                if whichfig == 'rsqrtfull':
                    if boxorplot == 'plot':
                        errorbarsplot(whichfig, data, void_name=voids_formal[draw], void_errname=voids_eformal[draw])  
                if whichfig == 'ratios':
                    print('inside fig')
                    theparams = 1  #two params
                    if boxorplot == 'box':
                        listofmeds = boxplotscatterone(whichfig, data, theparams=theparams)  
                    if boxorplot == 'plot':
                        errorbarsplot(whichfig, data, void_name=voids_formal[draw])  
                if whichfig == 'sigmas':
                    theparams = None  
                    if boxorplot == 'box':
                        listofmeds = boxplotscatterone(whichfig, data, theparams=theparams)  
                    if boxorplot == 'plot' and typeofsigma == None:
                        infdata = data[-1]/data 
                        errorbarsplot(whichfig, infdata, void_name=voids_formal[draw])  
                    if boxorplot == 'plot' and typeofsigma == 'free':
                        errorbarsplot(whichfig, data, void_name=voids_formal[draw])  
            
            if plotstyle == 'subpanel':
                whichfig = 'rsqrtfull'
                data = list(Fptr['%s_%s'%(whichfig,void_size)])  #returns a list 
                errorbarsplot(whichfig, data, void_name=voids_formal[draw], void_errname=voids_eformal[draw], styleplot=plotstyle)

                whichfigtwo = 'kvals' 
                figtwodata =list(Fptr['%s_%s'%(whichfigtwo,void_size)])  

                whichfig = 'rsqrttrue'
                data = list(Fptr['%s_%s'%(whichfig,void_size)])  #returns a list 
                errorbarsplot(whichfig, data, void_name=voids_formal[draw], void_errname=voids_eformal[draw], styleplot=plotstyle,
                              secondfig=whichfigtwo, seconddata=figtwodata)
                whichfig = 'ratios'
                data = list(Fptr['%s_%s'%(whichfig,void_size)])  #returns a list 
                errorbarsplot(whichfig, data, void_name=voids_formal[draw], void_errname=voids_eformal[draw], styleplot=plotstyle, 
                              secondfig=whichfigtwo, seconddata=figtwodata)


        # POWERSPECTRA: #make this a small definition too?? 
        if whichfig == 'kvals2D':
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
            x_offset = 0.08 #0.05
            y_offset = 4000 
            y_offset100 = 2500 

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

                for i, a in enumerate(target_snr):  #len of 7
                    if i == 0:
                        ax.plot(ps_xdata, medians_snr[i], color='r', alpha=1.0, linewidth=0.8)
                        continue
                    if i != len(medians_snr):
                        ax.plot(ps_xdata, medians_snr[i], color=colors[draw], alpha=1.0, linewidth=0.8)
                        ax.plot(ps_xdata, medians_noise[i], color='grey',alpha=1.0, linewidth=0.8, linestyle='dashed')
                        if a == 3:
                            plt.text(ps_xdata[-1] + x_offset, medians_noise[i][-1]+0.006, str(a), fontsize=12, ha='left', va='center') 
                            continue
                        if a == 5:
                            #plt.text(ps_xdata[-1] + x_offset, medians_noise[i][-1]-0.0002, str(a), fontsize=12, ha='left', va='center') 
                            plt.text(ps_xdata[-1] + x_offset, medians_noise[i][-1]-0.002, str(a), fontsize=12, ha='left', va='center') 
                            continue
                        if a == 10:
                            #plt.text(ps_xdata[-1] + x_offset, medians_noise[i][-1]-0.0003, str(a), fontsize=12, ha='left', va='center') 
                            plt.text(ps_xdata[-1] + x_offset, medians_noise[i][-1]-0.0004, str(a), fontsize=12, ha='left', va='center') 
                        else:
                            plt.text(ps_xdata[-1] + x_offset, medians_noise[i][-1], str(a), fontsize=12, ha='left', va='center') 

                    if i == len(target_snr[1:]): 
                        ax.plot(ps_xdata, medians_snr[i], color=colors[draw], alpha=1.0, linewidth=0.8, label='S+N') 
                        ax.plot(ps_xdata, medians_noise[i], color='grey',alpha=1.0, linewidth=0.8, linestyle='dashed', label='N') 
                        ax.plot(ps_xdata, medians_snr[0], color='r', alpha=1.0, linewidth=0.8, label='S') 
                        ax.fill_between(ps_xdata, the_lowers[0], the_uppers[0], color='r', alpha=0.2)
                    if i == 1: 
                        plt.text(ps_xdata[-1]+x_offset, medians_noise[i][-1]+y_offset100, 'SNR', fontsize=12, ha='left', va='center')  
                        #plt.text(ps_xdata[-1]+x_offset, medians_noise[i][-1]+y_offset, 'SNR', fontsize=12, ha='left', va='center')  

            # SAVING! 
            if k_plot == 'divkmax':  #for k/k_max
                labeling = r'$k/k_{\mathrm{max}}$'
            ax.set(xscale='log', yscale= 'log', xlim= (1e-2,1.7e0+x_offset), ylim= (1e-7,1e4+y_offset)) 
            ylabel=r'$P_{\mathrm{2D}}$'
            ax.set_xlabel(labeling, fontsize=12)
            ax.set_ylabel(ylabel, fontsize=12)
            #ax.legend(loc='upper left', fontsize=12, labelspacing=0.1, bbox_to_anchor=(0.01, 1.01))
            ax.legend(loc='upper left', fontsize=12, labelspacing=0.1, bbox_to_anchor=(0.01, 1.03))  #for 100um
            plt.savefig('sort_plots/boxplot15_powspec14_i1k_%s_%s.pdf'%(k_plot,void_size))
            print('saved spectra14!')
            plt.clf()
            Fptr.close()         

            # A QUICK CHECK TO SEE IF P2D VS K LOG-LOG POWERLAW IS INDEED 1. YES. 
            if 0:
                quickx = np.log10(ps_xdata)
                quicky = np.log10(nydata)
                rets = cfp.fit(linfunct, quickx[1:], quicky[-1][1:])
        Fptr.close() 

    if whichfig != 'kvals2D':
        plotfig(plotstyle) 



# plotting sigma rho calculated vs sigma rho brunt 
if 0:
    voids = ['50mu', '100mu', '12p5mu']
    colores = ['b', 'g', 'orange']
    for num, void_size in enumerate(voids):
        hfivefile = '/home/552/lv9897/sflaser/h5files/15sigma3D_%s.h5'%void_size
        hfivefile_B = '/home/552/lv9897/sflaser/h5files/15sigmaBrunt_%s.h5'%void_size
        Fptr = h5py.File(hfivefile,'r')
        Fptr_B = h5py.File(hfivefile_B,'r')
        ydata_3D = list(Fptr['sigma_3Ds_%s'%void_size])  
        ydata_Brunt = list(Fptr_B['sigma_brunts_%s'%void_size])  

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



# a very rough test with experimental data - make own file - DELETE
if 0:
    #fig,ax = plt.subplots(1,1)
    xdata = np.array(target_snr[:-1])

    teatime = ['t1', 't2']
    if 0:
        for tea in teatime:
            analyze = mpimg.imread('r60_%s_preshockreg.png'%tea)
            kvals, p2d, dk = main_snr_getdata.get_spectrumcf(analyze)
            p3d = 2 * kvals * p2d

            P2D = ((p2d[1:] * dk).sum()) *dk  #avoiding dk inside does not alter oom of sigmas
            P3D = ((p3d[1:] * dk).sum()) *dk

            sigmaNdk = np.sqrt(P2D)
            sigmaN = np.sqrt((p2d[1:] * dk).sum()) 

            sigmaRho_sqdk = sigmaNdk * (P3D/P2D)
            sigmaRho_sq = sigmaN * (P3D/P2D)
            
            sigmaRhodk = np.sqrt(sigmaRho_sqdk)
            sigmaRho = np.sqrt(sigmaRho_sq)

            plt.scatter(xdata, ydata)
    
    if 0:
        sigmasBs_dcc = [16.8, 16.2, 16.5] #uniform, 60, 120
        ydata = [sigmasBs_dcc[2]] * len(xdata)
        plt.scatter(xdata, ydata, c='k', alpha=0.5)

    if 0:
        plt.yscale='log'
        plt.xscale='log'
        thevoids = ['12p5mu', '50mu', '100mu']
        plt.savefig('quickbruntsnrtest_dccexp_%s.png'%thevoids[2])
        print('%s_saved!')
        plt.clf()


'''
# PLAN - how things got started: delete from here soon.
i. slice of positions of the cubes 
ii. plot 2d power spectrums of noise only and signal + noise

# for any box(position, proj, void) for each SNR 
a. get R^{1/2}_full values, where final is R^{1/2}_pure - call SNR function(k values turned off)

b. find the k values: AUTOMATE [DONE] 
    1. based on ii, for a possible range of k values for a specific SNR target, boxplot R^{1/2} vs that k range; iterations = 1000
    2. pick the k value that matches R^{1/2} with R^{1/2}_inf most closely BY HAND; store 
# pseudocode on automating
index = np.argmin(np.abs( r^1/2array - r1/2value ))
kcross = k_crossarray[index]
# these k values to recover rsqrt_true, these were picked by hand, could now check matching 
#ksnr = [6, 8, 18, 29, 35, 44, 52]  #right, left:[5, 7, 16, 25, 30, 38, 52], center: ksnr = [8, 8, 9, 19, 30, 36, 44, 52]

c. get R^{1/2}_true values - call SNR function(k values turned on)

# figures for any box(position, proj, void) or boxes combined 'superplot'
# plot or boxplot ...
1. R^{1/2}_full vs SNR
2. R^{1/2}_true vs SNR with R^{1/2}_pure
3a. R^{1/2}_true/R^{1/2}_full vs SNR
3b. Fit 3a. - Make a table with all the parameters
4?. k values vs SNR
'''
