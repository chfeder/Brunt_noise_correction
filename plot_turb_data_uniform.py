#!/usr/bin/env python

from globals_turb_data import *
from itertools import product
import argparse

# ===================================================================
def load_default_plotstyle():
    import matplotlib as mpl
    cfp.unload_plot_style()
    cfp.load_plot_style(fontsize=14)
    mpl.rcParams.update({
        'errorbar.capsize': 3,
        'lines.markersize': 5
    })

# ===================================================================
def plot_spectra(turb_data, flow_type: str, mach, noise_exp=1.0, res=256,
                res_resample=1.0, normalize=False, ax=None, save_panel=True):
    xlim = [1e-2, 2]
    ylim = [1e-14, 2e-2]
    if normalize:
        ylim = [1e-8, 3e4]
    kmax = turb_data['stats'][0].avg.k.max()
    x = turb_data['stats'][0].avg.k[1:-1]/kmax
    for isnr, SNR in enumerate(SNR_list):
        pos_correct = -0.2
        if SNR > 1e2: continue
        if SNR ==  3: pos_correct = +0.02
        if SNR ==  5: pos_correct = -0.32
        if SNR == 10: pos_correct = -0.50
        legend_xpos, legend_ypos = 0.015, 0.925
        legend_length = 1.35
        if SNR == 1: labels = ['N', 'S+N', r'$k_\mathrm{noise}/k_\mathrm{max}$']
        else: labels = [None, None, None]
        norm = 1.0
        if normalize:
            norm = turb_data['stats'][isnr].avg.PS_S[1:].sum()
        # draw the noise spectra
        y = turb_data['stats'][isnr].avg.PS_N[1:-1]/norm
        if SNR == 1: lf = cfp.legend_formatter().add(pos=(legend_xpos,legend_ypos-0.07), length=legend_length)
        else: lf = None
        cfp.plot(x=x, y=y, ax=ax, color='grey', linestyle=(0, (2, 2)), linewidth=1.0, label=labels[0], legend_formatter=lf)
        cfp.plot(x=1.1, y=10**(np.log10(y[-1])+pos_correct), ax=ax, text=str(SNR))
        # draw the noisy data spectra
        y = turb_data['stats'][isnr].avg.PS_SpN[1:-1]/norm
        if SNR == 1:
            lf = cfp.legend_formatter().add(pos=(legend_xpos,legend_ypos), length=legend_length)
        else:
            lf = None
        po = cfp.plot(x=x, y=y, ax=ax, color=flow_colors[flow_type], label=labels[1], legend_formatter=lf)
        # draw k_noise
        k_noise = turb_data['stats'][isnr].avg.k_noise/kmax
        k_noise_min = turb_data['stats'][isnr].p16.k_noise/kmax
        k_noise_max = np.min([turb_data['stats'][isnr].p84.k_noise/kmax, 1.0])
        if SNR > 1e-2 and k_noise_min >= 2.5e-2 and k_noise <= 0.93: # exclude drawing very high SNR cases, as k_noise is effectively 1
            y_cross = turb_data['stats'][isnr].avg.PS_S_at_k_noise/norm # intersection with signal spectrum
            po.plt().vlines(k_noise, ylim[0], y_cross, color=flow_colors[flow_type], linestyle=(0, (1, 1)))
            po.plt().fill_between([k_noise_min, k_noise_max], ylim[0], y_cross,
                                    color=flow_colors[flow_type], alpha=0.2, edgecolor='none')
            # fake line to get legend formatter to work
            if SNR == 1:
                lf = cfp.legend_formatter().add(pos=(legend_xpos+0.25,legend_ypos), length=legend_length)
            else:
                lf = None
            cfp.plot(x=0, y=1, ax=ax, color=flow_colors[flow_type], linestyle=(0, (1, 1)),
                        label=labels[2], legend_formatter=lf)
    # plot signal spectrum
    y = turb_data['stats'][0].avg.PS_S[1:-1]/norm
    yerr = [y-turb_data['stats'][0].p16.PS_S[1:-1]/norm, turb_data['stats'][0].p84.PS_S[1:-1]/norm-y]
    cfp.plot(x=x, y=y, ax=ax, yerr=yerr, shaded_err=True, color='red', label='S',
                legend_formatter=cfp.legend_formatter().add(pos=(legend_xpos,legend_ypos-0.14), length=legend_length))
    # generate final plot and save
    cfp.plot(x=0.886, y=0.92, ax=ax, text='SNR', normalised_coords=True)
    noise_str = ""
    if noise_exp != 1.0:
        noise_str = "_noise_exp_"+str(noise_exp)
    save = None
    if save_panel:
        save = figpath+'PS_N'+str(res)+flow_type+make_mach_str(mach)+noise_str+'.pdf'
    cfp.plot(xlabel=r"$k/k_\mathrm{max}$", ylabel=r"$P$", ax=ax, xlog=True, ylog=True,
             xlim=xlim, ylim=ylim, save=save)
# ===================================================================
def plot_k_noise_vs_SNR_for_flows(turb_data, flow_type, ax=None, save_panel=True):
    # ----------------------------------------------------
    def k_noise_func(SNR, SNR_crit=1.0):
        PL_exp = 2.0 / 3.0
        sharpness = 4.0
        k_noise = (1 / (1 + (SNR/SNR_crit)**(-PL_exp*sharpness)))**(1/sharpness)
        return np.log(k_noise)
    # ----------------------------------------------------
    # similar to figure 5 (k_noise vs SNR), but for fixed 50mu voids and resolution resampling factors 0.5, 1.0, 2.0
    xlim = [5e-2, 2e2]
    mach_list = []
    SNR_crit = []
    SNR_crit_err = []

    for iMach, mach in enumerate(mach_nums):
        kmax = turb_data[flow_type][mach]['stats'][0].avg.k.max()
        mach_list.append(mach)
        snr = []
        k_noise = []
        k_noise_err = []
        for isnr, SNR in enumerate(SNR_list):
            if SNR > 100 or SNR < 1e-1:
                continue
            x = SNR
            snr.append(x)
            y = turb_data[flow_type][mach]['stats'][isnr].avg.k_noise/kmax
            k_noise.append(y)
            yerr = [y-turb_data[flow_type][mach]['stats'][isnr].p16.k_noise/kmax,
                    turb_data[flow_type][mach]['stats'][isnr].p84.k_noise/kmax-y]
            k_noise_err.append(np.mean(yerr))
            if SNR == 1:
                label = r"$Mach="+str(mach)+r"$"
                lf = cfp.legend_formatter().add(pos=(0.64-(2-0)*0.05,0.45-iMach*0.085))
            else:
                label = None
                lf = None
            cfp.plot(x=x, y=y, yerr=np.array([yerr]).T, color=mach_colours[iMach], linestyle='none',
                        marker=flow_markers[flow_type], alpha=0.8, label=label, legend_formatter=lf, ax=ax)
        snr = np.array(snr)
        k_noise = np.array(k_noise)
        k_noise_err = np.array(k_noise_err)
        params = {'SNR_crit':[0.1,1.0,50]}
        fitdat = cfp.fit(k_noise_func, snr, np.log(k_noise), yerr=k_noise_err/k_noise,
                            params=params, n_random_draws=100)
        SNR_crit.append(fitdat.popt[0])
        SNR_crit_err.append(fitdat.perr[0])
        xfit = cfp.logspace(xlim[0], xlim[1], 200)
        yfit = np.exp(k_noise_func(xfit, *fitdat.popt)) # best-fit line

        lf = cfp.legend_formatter().add(pos=(0.01,0.92-0*0.075))

        cfp.plot(x=xfit, y=yfit, ax=ax, linewidth=1.0, color=mach_colours[iMach], linestyle=flow_linestyles[flow_type],
                    label=label, legend_formatter=lf)
    mach_list = np.array(mach_list)
    SNR_crit = np.array(SNR_crit)
    SNR_crit_err = np.array(SNR_crit_err)
    # generate final plot and save
    save = None
    if save_panel:
        save = figpath+"flow_type_"+str(flow_type)+'_k_noise_vs_SNR.pdf'
    cfp.plot(xlabel=r"$\mathrm{SNR}$", ylabel=r"$k_\mathrm{noise}/k_\mathrm{max}$", ax=ax, save=save,
                xlog=True, ylog=True, xlim=xlim, ylim=[2e-2,2e0])


# Creates a flow_types x machs grid of power spectra plots
def plot_spectra_grid(turb_data, flow_types=['s', 'm', 'c'], machs=[0.1, 0.5, 2, 10]):
    fig, axes = plt.subplots(len(machs), len(flow_types), figsize=(11, 3), sharex=True, sharey=True, squeeze=False)
    for icol, ft in enumerate(flow_types):
        for irow, mach in enumerate(machs):
            plot_spectra(turb_data[ft][mach], ft, mach, ax=axes[irow, icol], save_panel=False)
    fig.tight_layout()
    fig.savefig(figpath + 'PS_grid.pdf', bbox_inches='tight')
    plt.close(fig)

def plot_flows_grid(turb_data, flows=['s','c']):
    fig, axes = plt.subplots(len(flow_types), 1, figsize=(14, 9), sharex=True, sharey=True, squeeze=False)
    for irow, ft in enumerate(flows):
            plot_k_noise_vs_SNR_for_flows(turb_data, ft, ax=axes[0, irow], save_panel=False)
    fig.tight_layout()
    fig.savefig(figpath + 'Flows_grid.pdf', bbox_inches='tight')


# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    # Parse the -isim argument (if it exists).
    # This is used for Parallelising, by assigning each sim to a different core.
    parser = argparse.ArgumentParser()
    parser.add_argument("-isim", type=int, default=None, help="index of single sim to process (for job arrays)")
    args = parser.parse_args()

    load_default_plotstyle()

    sims = list(product(resolutions, flow_types, mach_nums))
    # Reduce the list of simulations to do to just the 1 per core if isim is used
    if args.isim is not None:
        sims = [sims[args.isim]]

    # load default data
    turb_data = {flow_type: {} for flow_type in flow_types}

    for res, flow_type, mach in sims:
        sim_name = make_model_name_params(res, flow_type, mach)
        dillfile = f"{project_path}processed_data/{sim_name}.dill"
        with open(dillfile, "rb") as f:
            turb_data[flow_type][mach] = dill.load(f)
        print("Read '"+dillfile+"'.", color='green')


    for mach in [0.1, 1, 10]:
        for ft in ['s', 'c']:
            print(ft, mach, [turb_data[ft][mach]['stats'][i].avg.k_noise for i in range(1,8)])


    # -------- FIGURES NOT TO BE USED IN REPORT --------

    """
    # Power spectrum plotted per sim
    for res, flow, mach in sims:
       plot_spectra(turb_data[flow][mach], flow, mach, res=res)


    # Machs plotted on k_noise vs SNR, plot for each flow
    for _, flow, _ in sims:
        plot_k_noise_vs_SNR_for_flows(turb_data, flow_type=flow)
    """


    # --------- FIGURES TO BE USED IN REPORT ------------

    # Grid of power spectra
    plot_spectra_grid(turb_data, flow_types, machs=[0.1, 0.5, 2, 10])

    # Grid of plots of x=SNR, y=knoise/kmax, across mach numbers
    plot_flows_grid(turb_data, flows=['c', 's'])

