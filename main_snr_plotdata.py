
from startlaser import * 
from collections import defaultdict
#from labelLines import *
# GET INFO FROM STORAGE 


# POSSIBLE FIT FUNCTIONS
def linfunct(x, a, b=1):  
    y = a * x + b
    return y

#def lininvfunct(x, a, b=1):
def lininvfunct(x, a, b):
    y = a * x**(-1) + b
    return y

def tanfunct_simple(x, first, second):  #used to be x, c, d, e; try this first when data is ready
    y = first * np.tanh(second*np.log10(x)) + (1-first)  #c should be negative, d positive; c+e=1
    return y

def sigmarho_fit(x,a):  
    #y = 1/ (lininvfunct(x, a, b=1.0) * tanfunct_simple(x, first=first, second=second)) 
    #y = (lininvfunct(x, a, b=1) * tanfunct_simple(x, first=first, second=second)) 
    pre_y = (lininvfunct(x, a, b=1) * tanfunct_simple(x, first=first, second=second)) 
    y = np.log(pre_y)
    return y 


# GET THE FITTING 
def fitthebox(whichfigure, eldata, secondparam):  
    '''
    GET the median, the 16th ('-' correction), and the 84th ('+' correction) of fit params a, b, c, d, etc
    this then allows us to report param +/- correction 
    '''
    param_a = []
    param_b = []
    param_c = []
    param_d = []
    param_e = []
    for j in range(numofboxes):  #avoid hard-coding numbers, but there are 15 boxes so...
        oneboxproj = []
        for k in range(len(xdata)):
            prexdata = eldata[k][j]
            oneboxproj = np.append(oneboxproj, prexdata)  #sigma_rho(SNR->inf)/sigma_rho(SNR) 
        
        if whichfigure == 'ratios': 
            print('inside the ratios')
            params = {"first":[-2,-1,-0.1], "second":[0.5,2.0,4]} 
            rets = cfp.fit(tanfunct_simple, xdata, oneboxproj, params=params) 
        if whichfigure == 'sigmas':
            print('inside the sigmas')
            params= {'a':[1e-10, 1, 1e10]}
            #secondparam = 1 #EDIT!
            rets = cfp.fit(sigmarho_fit, xdata, np.log(oneboxproj), params=params)#, yerr=)  #relative error 
            
        param_a = np.append(param_a, rets.popt[0])
        if secondparam == 1:
            param_b = np.append(param_b, rets.popt[1])
        if 0:
            param_c = np.append(param_c, rets.popt[2])
        if 0:
            param_d = np.append(param_d, rets.popt[3])
        if 0:
            param_e = np.append(param_e, rets.popt[4])

    med_a = np.percentile(param_a, 50) #same as np.median()
    p16_a = np.percentile(param_a, 16)
    p84_a = np.percentile(param_a, 84)
    if secondparam == 1:
        med_b = np.percentile(param_b, 50)
        p16_b = np.percentile(param_b, 16)
        p84_b = np.percentile(param_b, 84)
    if 0:
        med_c = np.percentile(param_c, 50)
        p16_c = np.percentile(param_c, 16)
        p84_c = np.percentile(param_c, 84)
    if 0:
        med_d = np.percentile(param_d, 50)
        p16_d = np.percentile(param_d, 16)
        p84_d = np.percentile(param_d, 84)
    if 0:
        med_e = np.percentile(param_e, 50)
        p16_e = np.percentile(param_e, 16)
        p84_e = np.percentile(param_e, 84)

    if secondparam == 0:
        med_b = p16_b = p84_b = 0
    if 1: 
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
        param_meds = [med_a, med_b]
        param_16s = [p16_a, p16_b]#, p16_c, p16_d, p16_e]  #report these as - for param
        param_84s = [p84_a, p84_b]#, p84_c, p84_d, p84_e]  #report these as + for param
    
    return param_meds, param_16s, param_84s


# TO PLOT THE FITTING 
def plotthefitting(thefigures, themedians, sixteenths, eightyfours):
    if thefigures == 'sigmas':
        the_yrecipemed = sigmarho_fit(the_xrecipe, *themedians)  
        the_yrecipelow = sigmarho_fit(the_xrecipe, *sixteenths)
        the_yrecipeupp = sigmarho_fit(the_xrecipe, *eightyfours)
    if thefigures == 'ratios':
        the_yrecipemed = tanfunct_simple(the_xrecipe, *themedians)  
        the_yrecipelow = tanfunct_simple(the_xrecipe, *sixteenths)
        the_yrecipeupp = tanfunct_simple(the_xrecipe, *eightyfours)
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
        the_data = ladata/ladata[-1]
        data_lista = the_data.tolist()
    if thefigure == 'ratios':
        data_lista = ladata[1:]
    else:
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

    if thefigure in {'sigmas'}:
        # call the fittings procedure definition
        meds, lows, upps = fitthebox(thefigure, data_lista, theparams)
        # call the fitting definition, which will provide the the_yrecipes below
        ln_ymed, ln_ylow, ln_yupp = plotthefitting(thefigure, meds, lows, upps)
        ymed = np.exp(ln_ymed)
        ylow = np.exp(ln_ylow)
        yupp = np.exp(ln_yupp)
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
        ax.set_ylabel('k noise')
        ax.set_yscale('log')
    if thefigure == 'sigmas':
        ax.set_yscale('log')
        ax.set_ylabel(r'$\sigma_{\rho/\rho_o,snr}$')# / \sigma_{\rho/\rho_o,inf}$') 
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




# MAIN: LOOP FOR THREE VOIDS!
numofboxes = 15
target_snr = [0.01, 0.1, 1, 3, 5, 10, 100, inf]  
subxdata = np.array(target_snr[2:-1])
the_subxrecipe = np.linspace(subxdata.min(), subxdata.max(), num=2500)  

whichfig = 'kvals2D'
if whichfig == 'ratios':
    xdata = np.array(target_snr[1:-1])
else:
    xdata = np.array(target_snr[:-1])
the_xrecipe = np.linspace(xdata.min(), xdata.max(), num=2500)  

k_plot = 'knoise'  #or divkmax or knoise
whichfigtwo = None#'rsqrttrue'
mainfigure = None#'sigmas'
plotstyle = 'singular'

#voids = ['50mu', '100mu', '12p5mu']
voids = ['100mu']#, '100mu', '12p5mu']
colors = ['g', 'b', 'orange']

if 1:
    for draw, void_size in enumerate(voids):
        print('void size is ',void_size)
        hfivefile = '/home/552/lv9897/sflaser/h5files/boxplot15_powspec14_i1k_SintNk75_%s.h5'%void_size
        Fptr = h5py.File(hfivefile,'r')

        # ALL OTHER DATA: BOXPLOTS + SCATTERS + FITS
        if whichfig != 'kvals2D':
            # these are parameters for the tanfunct_simple fit, initialized to zero for multiple plot purposes 
            first = 0
            second = 0
            data = list(Fptr['%s_%s'%(whichfig,void_size)])  #returns a list 
            #stop()
            if whichfigtwo == 'rsqrttrue':
                data_two = list(Fptr['%s_%s'%(whichfigtwo,void_size)])  
                listofmeds = boxplotscatter(whichfig, data, theparams=1, secondfig=whichfigtwo, seconddata=data_two) 

            # keep it singular paneled
            elif plotstyle == 'singular':
                listofmeds = boxplotscatterone(whichfig, data, theparams=1) #0 if first fig is sigma 
            else: 
                listofmeds = boxplotscatter(whichfig, data, theparams=1)

            if whichfig == 'kvals':
                for val, kvalues in enumerate(data):
                    data[val] = data[val]/data[-1] 

            if mainfigure == 'sigmas':
                data_three = list(Fptr['%s_%s'%(mainfigure,void_size)])  
                first = listofmeds[0]
                second = listofmeds[1]
                if plotstyle == 'singular':
                    boxplotscatterone(mainfigure, data_three, theparams=0)
                else: 
                    boxplotscatter(mainfigure, data_three, theparams=0)
            


        # POWERSPECTRA: #make this a small definition too?? 
        if whichfig == 'kvals2D':
            fig,ax = plt.subplots(1,1)
            
            whichfigtwo = 'thepowspec' 
            whichfigthree = 'thenoisepowspec'  
            whichfigfour = 'thesignalpowspec'
            ps_xdata = Fptr['%s_%s'%(whichfig,void_size)][()]  #returns a numpy array
            if k_plot == 'divkmax':  #for k/k_max
                for val, kvalues in enumerate(ps_xdata):
                    ps_xdata[val] = ps_xdata[val]/ps_xdata[-1] 
            if whichfigtwo:
                ydata = Fptr['%s_%s'%(whichfigtwo,void_size)][()]  
            if whichfigthree:
                nydata = Fptr['%s_%s'%(whichfigthree,void_size)][()]  
            if whichfigfour:
                sydata = Fptr['%s_%s'%(whichfigfour,void_size)][()]  

            for i, a in enumerate(target_snr[:-1]):  #len of 7
                #to display multiple iterations (choose how many), then also do ydata[i][ii]
                ax.plot(ps_xdata, ydata[i], color=colors[draw], alpha=1.0, linewidth=0.7)#, label=str(a))
                ax.plot(ps_xdata, nydata[i], color='grey',alpha=0.8, linewidth=0.7, linestyle='dashed')
                if i == 6:  #multiple iters also do ii==19
                    ax.plot(ps_xdata, ydata[i], color=colors[draw], alpha=1.0, linewidth=0.7, label='S+N') 
                    ax.plot(ps_xdata, nydata[i], color='grey',alpha=0.8, linewidth=0.7, linestyle='dashed', label='N') 
                    ax.plot(ps_xdata, sydata, color='r', alpha=1.0, linewidth=0.7, label='S') 
            # to label each line, could try:
            #pip install matplotlib-label-lines
            #add label = str(a) to plot, then
            #labelLines(plt.gca().get_lines(),zorder=2.5)

            
            if k_plot == 'divkmax':  #for k/k_max
                labeling = r'$k/k_{\mathrm{max}}$'
            if k_plot == 'knoise':  #for k/k_max
                labeling = r'$k_{\mathrm{noise}}$'
            ax.set(xscale='log', yscale= 'log', xlabel=labeling, ylabel=r'$P_{\mathrm{2D}}$') 
            ax.legend(loc='lower right')
            plt.savefig('sort_plots/boxplot15_powspec14_i1k_k75_%s_%s'%(k_plot,void_size))
            print('saved spectra14!')
            plt.clf()
            Fptr.close()         

            # A QUICK CHECK TO SEE IF P2D VS K LOG-LOG POWERLAW IS INDEED 1. YES. 
            if 0:
                quickx = np.log10(ps_xdata)
                quicky = np.log10(nydata)
                rets = cfp.fit(linfunct, quickx[1:], quicky[-1][1:])
        Fptr.close() 




# a very rough test with experimental data - make own file
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
    
    if 1:
        sigmasBs_dcc = [16.8, 16.2, 16.5] #uniform, 60, 120
        ydata = [sigmasBs_dcc[2]] * len(xdata)
        plt.scatter(xdata, ydata, c='k', alpha=0.5)

    plt.yscale='log'
    plt.xscale='log'
    thevoids = ['12p5mu', '50mu', '100mu']
    plt.savefig('quickbruntsnrtest_dccexp_%s.png'%thevoids[2])
    print('%s_saved!')
    plt.clf()

    # 'by eye'
    # SNR = 0.06
    # Sigma (40-50mu): 0.29
    # Sigma (100-120mu): 0.3
    






# NEXT: make a class? 
# PLAN: delete from here soon.
'''
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

'''
if 0:  #working on the cfp tool
y_err_up = the_ytenupp - the_ytenmed 
y_err_low = the_ytenmed - the_ytenlow 
y_err_arr = [y_err_low,y_err_up] 
cfp.plot(x=the_xten, y=the_ytenmed, yerr=y_err_arr, label='params median', shaded_err=True, save='testingthis.png', ylog=True)
'''