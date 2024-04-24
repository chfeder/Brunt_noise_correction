
from startlaser import * 
from collections import defaultdict
# GET INFO FROM STORAGE 


# POSSIBLE FIT FUNCTIONS
def linfunct(x, a, b):  
    y = a * x + b
    return y

def lininvfunct(x, a, b=1):
    y = a * x**(-1) + b
    return y

def tanfunct_simple(x, first, second):  #used to be x, c, d, e; try this first when data is ready
    y = first * np.tanh(second*np.log10(x)) + (1-first)  #c should be negative, d positive; c+e=1
    return y

def sigmarho_fit(x,a):  
    y = 1/ (lininvfunct(x, a, b=1.0) * tanfunct_simple(x, first=first, second=second)) 
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
    for j in range(15):  #avoid hard-coding numbers, but there are 15 boxes so...
        oneboxproj = []
        for k in range(len(xdata)):
            prexdata = eldata[k][j]
            if whichfigure == 'sigmas':
                oneboxproj = np.append(oneboxproj, prexdata)  #sigma_rho(SNR->inf)/sigma_rho(SNR) 
                #oneboxproj = np.append(oneboxproj, eldata[-1][j]/prexdata)  #sigma_rho(SNR->inf)/sigma_rho(SNR) 
            else: 
                oneboxproj = np.append(oneboxproj, prexdata) 
        
        if whichfigure == 'sigmas':
            params=None
            rets = cfp.fit(sigmarho_fit, xdata, oneboxproj, params=params) 
        else: 
            params = {"first":[-2,-1,-0.1], "second":[0.5,2.0,4]} 
            rets = cfp.fit(tanfunct_simple, xdata, oneboxproj, params=params) 
            
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
        #stop()
    
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
def boxplotscatter(thefigure, ladata, thesecond=0, secondfig=None, datasecond=None):
    fig,ax = plt.subplots(1,1)

    # funny order to make the boxplots with a plot overlayed to work
    plt.xscale('log')
    plt.xlim(np.array(xdata).min()/1.20, np.array(xdata).max()*1.20)

    # the boxplot!  WARNING: this sets only for 'sigmas' option 
    if thefigure == 'sigmas':
        the_data = ladata[-1]/ladata
        data_lista = the_data.tolist()
    else:
        data_lista = ladata
    boxplot = ax.boxplot(data_lista, showfliers=True, showmeans=False, labels=xdata, positions=xdata) #\
                        #, boxprops=dict(color="none"), capprops=dict(color="none"), medianprops=dict(color="none"))  

    # ...and hacking the box plot 
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
    plt.scatter(xdata, medians, color=colors[draw], marker='o', label='Median', alpha=0.5)  #reduce the marker size! 

    # for overlaying purposes, such as rhalf noise, rhalf denoised; may end up not needing/otherwise revise/edit
    if secondfig: 
        medians_two=[]
        # FOR 15 boxes, all target
        if 0:
            for k in range(7):
                medians_two = np.append(medians_two, np.median(datasecond[k]))
            plt.scatter(xdata, medians_two, color='orange', marker='o', alpha=0.5, label='other median')  #reduce the marker size! 
        # FOR 15 boxes, one target - signal only
        if 0:
            medians_two = np.array([np.median(datasecond)]*len(xdata))
            plt.scatter(xdata, medians_two, color='orange', marker='o', alpha=0.5, label='other median')  #reduce the marker size! 


    # call the fittings procedure definition
    meds, lows, upps = fitthebox(thefigure, data_lista, thesecond)
    # call the fitting definition, which will provide the the_yrecipes below
    if thefigure =='ratios':
        ymed, ylow, yupp = plotthefitting(thefigure, meds, lows, upps)

        plt.plot(the_xrecipe, ymed, color=colors[draw], linestyle='dashed')
        plt.plot(the_xrecipe, ylow, color=colors[draw], linestyle='dotted', alpha=0.7, label='lower') 
        plt.plot(the_xrecipe, yupp, color=colors[draw], linestyle='dotted', alpha=0.7, label='upper')  
        plt.fill_between(the_xrecipe, ylow, yupp, color=colors[draw], alpha=0.2)

    # SAVE    
    # ALSO AUTOMATE SOOOOOON.
    #plt.ylabel(r'$\sigma_{\rho/\rho_o ,Br}$') 
    #plt.ylabel(r'$R^{1/2}_{noisy}$')
    #plt.ylabel(r'$R^{1/2}_{denoised}$, $R^{1/2}_{pure signal}$')
    #plt.ylabel(r'$R^{1/2}_{noise}$, $R^{1/2}_{denoised}$')
    if thefigure == 'kvals':
        #plt.yscale('log')
        plt.ylabel('k noise')
    if thefigure == 'sigmas':
        plt.yscale('log')
        plt.ylabel(r'$\sigma_{\rho/\rho_o,inf} / \sigma_{\rho/\rho_o,snr}$') 
    if thefigure == 'ratios':
        plt.ylabel(r'$R^{1/2}_{denoised}/R^{1/2}_{noisy}$')
    plt.xlabel('SNR')

    # really saving now
    plt.savefig('sort_plots/boxplot15_%s_i1kfix_%s'%(thefigure,void_size))
    print('%s_saved!'%thefigure)
    plt.clf()
    Fptr.close() 
    return meds




# MAIN: LOOP FOR THREE VOIDS!
target_snr = [0.01, 0.1, 1, 3, 5, 10, 100, inf]  
xdata = np.array(target_snr[:-1])
the_xrecipe = np.linspace(xdata.min(), xdata.max(), num=2500)  
voids = ['50mu', '100mu', '12p5mu']
colors = ['g', 'b', 'orange']
for draw, void_size in enumerate(voids):
    print('void size is ',void_size)
    whichfig = 'kvals'
    whichfigtwo = None 
    mainfigure = None #'sigmas'


    # ALL OTHER DATA: BOXPLOTS + SCATTERS + FITS
    if whichfig != 'kvals2D':
        hfivefile = '/home/552/lv9897/sflaser/h5files/boxplot15_i1kforwardk51ratiofix_%s.h5'%void_size
        Fptr = h5py.File(hfivefile,'r')
        # these are parameters for the tanfunct_simple fit, initialized to zero for multiple plot purposes 
        first = 0
        second = 0

        data = list(Fptr['%s_%s'%(whichfig,void_size)])  #returns a list 
        if whichfigtwo:
            data_two = list(Fptr['%s_%s'%(mainfigure,void_size)])  #returns a list 

        if whichfig == 'kvals':
            for val, kvalues in enumerate(data):
                data[val] = data[val]/data[-1] 
        boxplotscatter(whichfig, data, thesecond=1)  #see if it still works.
        if mainfigure == 'sigmas':
            first = listofmeds[0]
            second = listofmeds[1]
            boxplotscatter(mainfigure, data_two, thesecond=0)


    # POWERSPECTRA: #make this a small definition too?? 
    if whichfig == 'kvals2D':
        hfivefile = '/home/552/lv9897/sflaser/h5files/powspec14_i1kSNR_noisesignalwnoise_%s.h5'%void_size
        Fptr = h5py.File(hfivefile,'r')
        fig,ax = plt.subplots(1,1)
        
        whichfigtwo = 'thepowspec'
        whichfigthree = 'thenoisepowspec'
        ps_xdata = Fptr['%s_%s'%(whichfig,void_size)][()]  #returns a numpy array
        for val, kvalues in enumerate(ps_xdata):
            ps_xdata[val] = ps_xdata[val]/ps_xdata[-1] 
        if whichfigtwo:
            ydata = Fptr['%s_%s'%(whichfigtwo,void_size)][()]  
        if whichfigthree:
            nydata = Fptr['%s_%s'%(whichfigthree,void_size)][()]  

        for i in range(len(target_snr[:-1])):
            ax.plot(ps_xdata, ydata[i], color='blue', alpha=0.5, label='powspec') 
            ax.plot(ps_xdata, nydata[i], color='grey', label='powspec') 
        ax.set(xscale= 'log', yscale= 'log', xlabel='k values', ylabel=r'$P_{2D}$')#, xlim=(5,51)) #, ylim=(10e-1,10e11))#, xlim=(1e0,1e2))  #xscale= 'log',
        #ax.axvline(x=51)
        plt.savefig('sort_plots/powspec14_i1kSNR_divkmax_noisesignalwnoise_%s'%(void_size))
        print('saved spectra14!')
        plt.clf()
        Fptr.close()         




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
