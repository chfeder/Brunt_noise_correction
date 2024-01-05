
from startlaser import * 
from collections import defaultdict


# BUILDS THE BOX
def get_box(voids, box, projecting, proj_axis):
    the_ds = yt.load('/scratch/ek9/ccf100/nif/turb_foam_gamma_5_3_3d_1024_%s_bubbles/NIF_hdf5_plt_cnt_0375'%(voids))
    the_level = 0
    delta_a = 0.02  #(1.2/3)/2
    edge_mtop = 0.3675 
    edge_mbot = 0.3275

    if box == 'c_box':
        the_ledge = the_ds.arr([-delta_a, edge_mbot, -delta_a],'code_length')
        the_redge = the_ds.arr([delta_a, edge_mtop, delta_a],'code_length')

    if box == 'r_box':
        the_ledge = the_ds.arr([delta_a, edge_mbot, -delta_a],'code_length')
        the_redge = the_ds.arr([delta_a+2*delta_a, edge_mtop, delta_a],'code_length')
    if box == 'f_box':
        the_ledge = the_ds.arr([-delta_a, edge_mbot, delta_a],'code_length')
        the_redge = the_ds.arr([delta_a, edge_mtop, delta_a+2*delta_a],'code_length')

    if box == 'l_box':
        the_ledge = the_ds.arr([-delta_a-2*delta_a, edge_mbot, -delta_a],'code_length')
        the_redge = the_ds.arr([-delta_a, edge_mtop, delta_a],'code_length')
    if box == 'b_box':
        the_ledge = the_ds.arr([-delta_a, edge_mbot, -delta_a-2*delta_a],'code_length')
        the_redge = the_ds.arr([delta_a, edge_mtop, -delta_a],'code_length')

    the_center =0.5*the_ledge + 0.5*the_redge
    mindx = the_ds.index.get_smallest_dx()
    size = ((the_redge - the_ledge)/mindx).astype('int') 

    section = extract_s1.extract_section_1(the_ds,the_level,the_ledge,size)
    cg_first, cg_second, _ = project.Project(section, proj_axis, projecting, field_main, field_two)

    return cg_first, cg_second


def noise_iterations(iterations, cg_main, cg_sec, lowpassfilter, target, noise, ksnr=None):#, field_main='density'):
    collect = []
    for many in range(iterations):  

        # ADD NOISE TO THE IMAGE
        if noise == True:
            cg_wnoise, _ = noiseimage.add_noise(cg_main, target)
        if noise == False:
            cg_wnoise = cg_main
        
        # OBTAIN SIGMA N^2
        sigma_sq = calculate.Calculate(field_main, cg_wnoise, cg_sec) 

        # OBTAIN P2D AND P3D           
        if lowpassfilter == 'no':
            four_amps, powerspec2D, powerspec3D, infxyz, ks = powerspec.powerspec(cg_wnoise, dimension, target)    
        if lowpassfilter == 'yes':
            ks, powerspec2D, powerspec3D, infxyz = powerspec.powerspec_lowpass(cg_wnoise, the_k=ksnr, lowpass_filter=lowpassfilter)
        
        # OBTAIN SIGMA RHO OR ...
        sigma_rho, sigma_rho2d, brunts_rsq, _ = parsevals.parsevals_th(dimension, powerspec2D, powerspec3D, infxyz, ks, sigma_sq)
        #collect[target].append(brunts_rsq)   # what are you storing.
        collect = np.append(collect, brunts_rsq)

        if noise == False:
            break

    # OBTAIN MEAN VALUE 
    #the_rs = [l for l in collect.values()]
    #emptylist = []
    #the_rsall = np.append(emptylist,the_rs[val])

    the_rsallmed = np.median(collect)
    the_rsalllogmed = np.log10(the_rsallmed)
    the_rsalllog = np.log10(collect)
    the_rsallmedlog = np.median(the_rsalllog) 
    the_rsalllogerr = np.std(the_rsalllog)
    #stop() #logmed ~ medlog, check [ ]
    return the_rsallmed



# SETUP  
field_main = 'density'
field_two = 'cell_volume'
polytropic_index = 1 # +/-0.1, for 50mu voids

projection = True 
dimension = 1
noise = True
iterations = 5#1000
noiseat = 'image'  
lowpassfilter = 'no'  # yes:low k values, no:all k values

# SIGNAL TO NOISE RATIOs & COLORS: 
# for pure signal, simply do target_snr[-1], for all others target_snr[:-1]
target_snr = [0.01, 0.1, 1, 3, 5, 10, 100, inf]
target_colors = ['orange','y','magenta','blue','green','cyan','red'] 
# k values to recover rsqrt_true, these were picked by hand, to be automated: save to the bottom when it works or delete    
#ksnr = [6, 8, 18, 29, 35, 44, 52]  #right, left:[5, 7, 16, 25, 30, 38, 52], center: ksnr = [8, 8, 9, 19, 30, 36, 44, 52]


# ALL THE OPTIONS
# 3 different void sizes, but do for 1 only: 5 boxes x 3 projections: 15 cases
void_size = '50mu'

# 3 different box locations for the x los: left, center, right. and 2 more for the z los: front & back 
if 1:
    boxes = ['l_box', 'c_box', 'r_box', 'b_box', 'f_box'] 
if 0: #debug
    boxes = ['c_box'] 
# 3 different projections
if 1:  
    projection_axis_list = [0, 1, 2]  
if 0:#debug
    projection_axis_list = [0]



# START HERE
collect = defaultdict(list)
rsqrt_full =[]
rsqrt_pure =[]
rsqrt_true =[]
rsqrt_ratio =[]
kcross = []
k_array = np.linspace(1,75, 75)  #start, stop, num 

# three projection axes
for eje, los in enumerate(projection_axis_list):
    # five box locations
    for num, boxing in enumerate(boxes):  
        cg_main, cg_vol = get_box(void_size, boxing, projection, los)

        # seven SNR
        for val,target in enumerate(target_snr):
            lowpassfilter = 'no'  
            rsqrt_med = noise_iterations(iterations, cg_main, cg_vol, lowpassfilter, target, noise)
            if target != inf:
                rsqrt_full = np.append(rsqrt_full, rsqrt_med)
                #collect[target].append(rsqrt_med)  #box-plotting works 
            if target == inf:
                rsqrt_pure = np.append(rsqrt_pure, rsqrt_med)
     
        if 1:
            #print('rsqrt_full ',rsqrt_full)
            for val,target in enumerate(target_snr):
                rsqrt_k = []
                lowpassfilter = 'yes'  
                if target != inf:
                    for kval in k_array: 
                        rsqrt_med_k = noise_iterations(iterations, cg_main, cg_vol, lowpassfilter, target, noise, kval)
                        rsqrt_k = np.append(rsqrt_k, rsqrt_med_k)
                    index = np.argmin(np.abs(rsqrt_k - rsqrt_med))  #rsqrt_med should be the pure signal value here 
                    kcross = k_array[index] 
                    #print('kcross, box, proj ',kcross, num, los)
                    rsqrt_true = noise_iterations(iterations, cg_main, cg_vol, lowpassfilter, target, noise, kcross)
                    #print('true ',rsqrt_true)
                    rsqrt_ratio = rsqrt_true/rsqrt_full[val]
                    #print('full ',rsqrt_full[val])
                    #print('ratio ',rsqrt_ratio)

                    #rsqrt_ratio_log = np.log(rsqrt_ratio)
                    #target_log = np.log10(target)
                    collect[target].append(rsqrt_ratio) 


# POSSIBLE FIT FUNCTIONS
def cfunct(x, a, b, c, d):
    y = a * np.tanh((c*x)+d) + b
    return y

def dfunct(x, a, b):
    y = 1/(1 + a*np.exp(-x+b)) 
    return y

# BOXPLOTS WITH PLOTTINGS
fig,ax = plt.subplots(1,1)
if 1:
    xdata = target_snr[:-1] 
    xdatalog = np.log10(target_snr[:-1]) 

    labels, data = collect.keys(), collect.values()  
    boxplot = ax.boxplot(data, showfliers=True, showmeans=False)  #add some labels
    plt.close()

    medians = [median.get_ydata()[0] for median in boxplot['medians']]
    ydatalog = np.log10(medians)

    #plt.scatter(xdata, medians, color='red', marker='o', label='Median')   #this shrinks the plot
    if 0: #for boxplot
        plt.scatter(range(1, len(labels) + 1), medians, color='red', marker='o', label='Median')
    if 1: #for plot
        plt.scatter(xdata, medians, color='red', marker='o', label='Median')

    if 1:  #add the fit! 
        rets = cfp.fit(cfunct, xdatalog, ydatalog)
        the_xrecipe = np.linspace(xdatalog.min(),xdatalog.max(),num=500)  #more numbers (500) for the x fitting!;
                                                                                    #len(xdatalog) for a raw fit WITH boxplot
        the_xten = 10**the_xrecipe
        the_yrecipe = cfunct(the_xrecipe, *rets.popt)
        the_yten = 10**the_yrecipe
        plt.plot(the_xten, the_yten, 'g--')  #this yields a squeezed boxplot, but is a smooth plot.
        #plt.plot(range(1, len(labels) + 1), the_yten, 'g--')  #this gets it right, but the fit is boxy.
        print('param cf vals ',rets.popt) 

    plt.xlabel('SNR')
    plt.xscale('log')
    plt.ylabel(r'$R^{1/2}_{true}/R^{1/2}_{full}$')
    plt.yscale('log')
    #plt.xticks(range(1, len(labels) + 1), labels)
    #plt.savefig('testing_megaboxplot_15')
    plt.savefig('originalplot_15')#%(box,target,proj_axis,void_size))
    #plt.savefig('rsqrt_full_vs_SNR_%s_%d_%s'%(box,proj_axis,void_size))  #rsqrt_full, rsqrt_true, rhalf_ratio
    print('saved!')
    plt.clf()


# after the boxplots
if 0:  #tanh
    rets = cfp.fit(cfunct, xdatalog, ydatalog)
    the_xrecipe = np.linspace(xdatalog.min(),xdatalog.max(),num=500)  # more numbers for the x fitting! ,was len(xdatalog)
    the_xten = 10**the_xrecipe
    the_yrecipe = cfunct(the_xrecipe, *rets.popt)
    the_yten = 10**the_yrecipe
    plt.plot(the_xten, the_yten, 'g--')
    print('param cf vals ',rets.popt) 
    plt.scatter(xdata, ydata)  
    #plt.plot([1e-2,100],[rhalf_true,rhalf_true],'--',c='grey')
    #ax.set(ylim=[2e-1,1.1e0])  #ax.set(xlim=[1e-3,2e2], ylim=[1e-1,4e-1])





# SOME THINGS WE'D LIKE TO KNOW, SEE, ETC
# FOR PLOTS 
if 0:
    fig,ax = plt.subplots(1,1)
    fig.subplots_adjust(wspace=0, hspace=0)
# PLOT POWERSPECTRUMS
if 0:  
    k1, spec1, _, _ = powerspec.powerspec_lowpass(cg_wnoise)
    #k2, spec2, _, _ = powerspec.powerspectrum_ii(cg_first)
    k3, spec3, _, _ = powerspec.powerspec_lowpass(noiseonly, justnoise=noise)
    ax.plot(k1, spec1, c=target_colors[val], alpha=0.5)
    #ax.plot(k2, spec2, c='b')
    ax.plot(k3, spec3, c='k', alpha=0.5)
    ax.set(yscale='log')#, xscale='log')
# SAVE THE POWERSPECTRUM PLOTS 
if 0:  
    plt.ylim(1e1,1e11)
    #plt.xlim(40,75)
    fig.savefig('pd2_snr_sn_wn_logy_%s_%d.png'%(box,proj_axis))
# IMAGE THE SIGNAL, NOISE, SIGNAL PLUS NOISE
if 0:  
    cfp.plot_map(cg_first.v,save='signal_%s'%box)#, vmin=20, vmax=49)
    #cfp.plot_map(cg_wnoise,save='signal_noise.png')
    #cfp.plot_map(noiseonly,save='noise.png')
    #stop()


# PURPOSE 
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
c. get R^{1/2}_true values - call SNR function(k values turned on)

# figures for any box(position, proj, void) or boxes combined 'superplot'
# plot or boxplot ...
1. R^{1/2}_full vs SNR
2. R^{1/2}_true vs SNR with R^{1/2}_pure
3a. R^{1/2}_true/R^{1/2}_full vs SNR
3b. Fit 3a. - Make a table with all the parameters
4?. k values vs SNR
'''