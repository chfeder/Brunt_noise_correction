#!/usr/bin/env python

from globals_turb_data import *
from itertools import product
from datetime import datetime
import math
import argparse

# ===================================================================
# ------------------------ HELPER FUNCTIONS --------------------------
def load_default_plotstyle():
    import matplotlib as mpl
    cfp.unload_plot_style()
    cfp.load_plot_style(fontsize=14)
    mpl.rcParams.update({
        'errorbar.capsize': 3,
        'lines.markersize': 5
    })

# ============================= FIGURE 1 ========================================
def plot_z_proj_grid(turb_data, res, isnr=0, flow_types=['s','m','c'], machs=[0.1, 0.5, 2, 10], timestep=50,
                        grid_size=(9, 12), log=True, label_info=""):
    # j here is the number which refers to the ret object, since they follow a pattern of
    # (0,x) -> (0,y) -> (0,z) -> (1,x) -> ...
    # for (timestep, projection LOS).

    cfp.load_plot_style(fontsize=20)

    midpt = math.floor(len(machs)/2)

    lower_machs = machs[0:midpt]
    upper_machs = machs[midpt:len(machs)]

    # set the real min and max as log boundaries
    j = (timestep * 3) + 2

    if log:

        all_z_projs = [turb_data[ft][mach]['data'][isnr][j].noisy_data[0] for ft in flow_types for mach in lower_machs]
        lower_vmin = min(m.min() for m in all_z_projs)
        lower_vmax = max(m.max() for m in all_z_projs)
        print(f"lower_vmin: {lower_vmin}")
        print(f"lower_vmax: {lower_vmax}")

        all_z_projs = [turb_data[ft][mach]['data'][isnr][j].noisy_data[0] for ft in flow_types for mach in upper_machs]
        upper_vmin = min(m.min() for m in all_z_projs)
        upper_vmax = max(m.max() for m in all_z_projs)
        print(f"upper_vmin: {upper_vmin}")
        print(f"upper_vmax: {upper_vmax}")
    else:
        #original paper's: vmin, vmax = 5e-3, 1e-2
        lower_vmin, lower_vmax = 5e-3, 1e-2
        upper_vmin, upper_vmax = 5e-3, 1e-2


    cmap = 'magma'
    SNR = SNR_list_uniform[isnr]


    fig, axes = plt.subplots(len(machs), len(flow_types), squeeze=False, figsize=grid_size, gridspec_kw={'hspace': 0, 'wspace': 0})
    fig.set_constrained_layout(False)
    for irow, mach in enumerate(machs):
        for icol, ft in enumerate(flow_types):
            ax = axes[irow, icol]

            # check SNR match
            dat = turb_data[ft][mach]['data'][isnr][j]
            if dat.SNR[0] != SNR:
                print("SNR DOESN'T MATCH")

            if mach in lower_machs:
                po_lower = cfp.plot_map(dat.noisy_data[0], ax=ax, log=log, axes_format=["", ""], axes_ticks=[[],[]],
                                    vmin=lower_vmin, vmax=lower_vmax, cmap=cmap, colorbar=False)
                ax.set_box_aspect(1)
            else:
                po_upper = cfp.plot_map(dat.noisy_data[0], ax=ax, log=log, axes_format=["", ""], axes_ticks=[[],[]],
                                    vmin=upper_vmin, vmax=upper_vmax, cmap=cmap, colorbar=False)
                ax.set_box_aspect(1)

            ax.text(0.5, 0.88, r"$\mathcal{M} = "+str(mach)+r"$", ha='center', va='center',
                        transform=ax.transAxes, bbox=dict(boxstyle="round", fc="white", alpha=0.9))
            if irow == 0:

                ax.set_title(flow_labels[ft])

    #important: this solves the alignment problem, eliminating all the gaps
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0, hspace=0)

    save = figpath+f'main/fig_1_z_proj_grid_'+cfp.eform(SNR,prec=1)+f'_N{res}_{label_info}.pdf'
    fig.savefig(save)
    plt.close(fig)
    print("Wrote '"+save+"'.", color='magenta')


    #plot colorbar
    cfp.plot_colorbar(vmin=upper_vmin, vmax=upper_vmax, log=log, cmap=cmap, label=r'Column density $(\mathrm{g\,cm^{-2}})$',
                       panels=3, orientation="vertical", save=figpath+f'main/fig_1_cbar_above_zero_machs.pdf')
    cfp.plot_colorbar(vmin=lower_vmax, vmax=lower_vmax, log=log, cmap=cmap, label=r'Column density $(\mathrm{g\,cm^{-2}})$',
                       panels=3, orientation="vertical", save=figpath+f'main/fig_1_cbar_below_zero_machs.pdf')


    load_default_plotstyle()


# ============================= FIGURE 2 ====================================
# Creates a flow_types x machs grid of power spectra plots
def plot_spectra_grid(turb_data, res, flow_types=['s', 'm', 'c'], machs=[0.1, 0.5, 2, 10], grid_size=(12, 16)):
    fig, axes = plt.subplots(len(machs), len(flow_types), figsize=grid_size, sharex=True, sharey=True, squeeze=False)
    for icol, ft in enumerate(flow_types):
        for irow, mach in enumerate(machs):
            plot_spectra(turb_data[ft][mach], ft, mach, ax=axes[irow, icol], save_panel=False)

            # Add column and row titles
            if irow == 0:
                axes[irow, icol].set_title(flow_labels[ft])
            if icol == 0:
                axes[irow, icol].set_ylabel(r"$\mathcal{M} = " + str(mach) + "$")
    fig.tight_layout()
    fig.savefig(figpath + f'main/fig_2_PS_grid_N{res}.pdf', bbox_inches='tight')
    plt.close(fig)

def plot_spectra(turb_data, flow_type: str, mach, noise_exp=1.0, res=256,
                res_resample=1.0, normalize=False, ax=None, save_panel=True):
    xlim = [1e-2, 2]
    ylim = [1e-13, 1e4]
    if normalize:
        ylim = [1e-8, 3e4]
    kmax = turb_data['stats'][0].avg.k.max()
    x = turb_data['stats'][0].avg.k[1:-1]/kmax
    for isnr, SNR in enumerate(SNR_list_uniform):
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
        save = figpath+'extra/fig_2_PS_N'+str(res)+flow_type+make_mach_str(mach)+noise_str+'.pdf'
    cfp.plot(xlabel=r"$k/k_\mathrm{max}$", ylabel=r"$P$", ax=ax, xlog=True, ylog=True,
             xlim=xlim, ylim=ylim, save=save)


# =============================FIGURE 3 =======================================
def plot_k_noise_vs_SNR_for_flows(turb_data, flow_type, ax=None, save_panel=True):
    # ----------------------------------------------------
    def k_noise_func(SNR, SNR_crit=1.0, PL_exp=2/3):
        #PL_exp = 2.0 / 3.0
        sharpness = 4.0
        k_noise = (1 / (1 + (SNR/SNR_crit)**(-PL_exp*sharpness)))**(1/sharpness)
        return np.log(k_noise)
    # ----------------------------------------------------
    # similar to figure 5 (k_noise vs SNR), but for fixed 50mu voids and resolution resampling factors 0.5, 1.0, 2.0
    xlim = [5e-2, 2e4]
    mach_list = []
    SNR_crit = []
    SNR_crit_err = []
    PL_exp_list = []
    PL_exp_err = []

    for iMach, mach in enumerate(mach_nums):
        kmax = turb_data[flow_type][mach]['stats'][0].avg.k.max() #this should be context for the uniform experiemtn

        print(f"kmax is {kmax}")

        mach_list.append(mach)
        snr = []
        k_noise = []
        k_noise_err = []
        for isnr, SNR in enumerate(SNR_list_uniform):
            x = SNR
            snr.append(x)
            y = turb_data[flow_type][mach]['stats'][isnr].avg.k_noise/kmax
            k_noise.append(y)
            yerr = [y-turb_data[flow_type][mach]['stats'][isnr].p16.k_noise/kmax,
                    turb_data[flow_type][mach]['stats'][isnr].p84.k_noise/kmax-y]
            k_noise_err.append(np.mean(yerr))
            if SNR == 1:
                label = r"$Mach="+str(mach)+r"$"
                lf = cfp.legend_formatter().add(pos=(0.64-(2-0)*0.05,0.55-iMach*0.085))
            else:
                label = None
                lf = None
            cfp.plot(x=x, y=y, yerr=np.array([yerr]).T, color=mach_colours[iMach], linestyle='none',
                        marker=flow_markers[flow_type], alpha=0.8, label=label, legend_formatter=lf, ax=ax)
        snr = np.array(snr)
        k_noise = np.array(k_noise)
        k_noise_err = np.array(k_noise_err)

        # free variables to fit the function too
        # format: [lower bound, starting guess, upper bound]
        params = {
            'SNR_crit': [0, 1.0, 1000], #original: [0.1, 1.0, 50]
            'PL_exp':   [0, 1.0, 1000]  #originally fixed to 2/3
        }
        fitdat = cfp.fit(k_noise_func, snr, np.log(k_noise), yerr=k_noise_err/k_noise,
                            params=params, n_random_draws=100)

        SNR_crit.append(fitdat.popt[0])
        SNR_crit_err.append(fitdat.perr[0])
        xfit = cfp.logspace(xlim[0], xlim[1], 200)
        yfit = np.exp(k_noise_func(xfit, *fitdat.popt)) # best-fit line

        # I removed label=label and legend_formatter=lf here
        cfp.plot(x=xfit, y=yfit, ax=ax, linewidth=1.0, color=mach_colours[iMach], linestyle=flow_linestyles[flow_type])
    mach_list = np.array(mach_list)
    SNR_crit = np.array(SNR_crit)
    SNR_crit_err = np.array(SNR_crit_err)
    # generate final plot and save
    save = None
    if save_panel:
        save = figpath+"extra/flow_type_"+str(flow_type)+'_k_noise_vs_SNR.pdf'
    cfp.plot(xlabel=r"$\mathrm{SNR}$", ylabel=r"$k_\mathrm{noise}/k_\mathrm{max}$", ax=ax, save=save,
                xlog=True, ylog=True, xlim=xlim, ylim=[2e-2,2e0])


def plot_flows_grid(turb_data, res, flow_types=['s','c'], grid_size=(12, 16)):
    fig, axes = plt.subplots(1, len(flow_types), figsize=grid_size, sharex=True, sharey=True, squeeze=False)
    for icol, ft in enumerate(flow_types):
            axes[0, icol].set_title(flow_labels[ft])
            plot_k_noise_vs_SNR_for_flows(turb_data, ft, ax=axes[0, icol], save_panel=False)
    fig.tight_layout()
    fig.savefig(figpath + f'main/fig_3_flows_grid_N{res}.pdf', bbox_inches='tight')



# =========================== FIGURE 4 ==============================================
# somehwat of a copy of figure 3's, except this function doesn't plot, it only fits the data.
# used to create the subpanels.
def fit_k_noise_vs_SNR(turb_data, flow_type):
    # ----------------------------------------------------
    def k_noise_func(SNR, SNR_crit=1.0, PL_exp=2/3):
        sharpness = 4.0
        k_noise = (1 / (1 + (SNR/SNR_crit)**(-PL_exp*sharpness)))**(1/sharpness)
        return np.log(k_noise)
    # ----------------------------------------------------
    mach_list = []
    SNR_crit = []
    SNR_crit_err = []
    PL_exp_list = []
    PL_exp_err = []

    for iMach, mach in enumerate(mach_nums):
        kmax = turb_data[flow_type][mach]['stats'][0].avg.k.max()
        mach_list.append(mach)
        snr, k_noise, k_noise_err = [], [], []
        for isnr, SNR in enumerate(SNR_list_uniform):
            snr.append(SNR)
            y = turb_data[flow_type][mach]['stats'][isnr].avg.k_noise/kmax
            k_noise.append(y)
            yerr = [y-turb_data[flow_type][mach]['stats'][isnr].p16.k_noise/kmax,
                    turb_data[flow_type][mach]['stats'][isnr].p84.k_noise/kmax-y]
            k_noise_err.append(np.mean(yerr))
        snr = np.array(snr); k_noise = np.array(k_noise); k_noise_err = np.array(k_noise_err)

        # free variables to fit the function too
        # format: [lower bound, starting guess, upper bound]
        params = {
            'SNR_crit': [0, 1.0, 1000], #original: [0.1, 1.0, 50]
            'PL_exp':   [0, 1.0, 1000]  #originally fixed to 2/3
        }
        fitdat = cfp.fit(k_noise_func, snr, np.log(k_noise), yerr=k_noise_err/k_noise,
                         params=params, n_random_draws=100)

        SNR_crit.append(fitdat.popt[0]); SNR_crit_err.append(fitdat.perr[0])
        PL_exp_list.append(fitdat.popt[1]); PL_exp_err.append(fitdat.perr[1])

    return (np.array(mach_list), np.array(SNR_crit), np.array(SNR_crit_err),
            np.array(PL_exp_list), np.array(PL_exp_err))

def plot_fit_params_grid(turb_data, flow_types=['s', 'm', 'c']):
    """
    2-panel graph as follows:


    SNR_crit        |   1   |
    PL_exp = \alpha |   2   |
    ---------------------------------
                       Mach

    each panel has the param as the y-axis, and all have the mach as x-axis

    """
    fig, axes = plt.subplots(2, 1, figsize=(6, 8), squeeze=False)


    #def PL(x, a, b):
    #    return a*x**b

    def PL_log(x, log_a, b):
        return log_a + b * np.log(x)

    # Text placement constants
    y_of_highest_equality = 0.2 #equality, as in, the labels like "SNR = (3 +/- 4) x ..."
    x_of_upper_row_flow_type = 0.6


    print("Entering fig 4 loop")
    for i_flow_type, flow_type in enumerate(flow_types):
        mach_list, SNR_crit, SNR_crit_err, PL_exp_list, PL_exp_err = fit_k_noise_vs_SNR(turb_data, flow_type)

        icol = 0

        # --- row 0: SNR_crit vs Mach ---
        ax = axes[0, icol]
        cfp.plot(x=mach_list, y=SNR_crit, yerr=SNR_crit_err.T, ax=ax,
                    marker=flow_markers[flow_type], linestyle='none', color=flow_colors[flow_type],
                    label=flow_labels[flow_type],
                    legend_formatter=cfp.legend_formatter().add(pos=(x_of_upper_row_flow_type,0.90-i_flow_type*0.09))
                    )
        # Log fit code
        yerr_mean = np.abs(SNR_crit_err.T).mean(axis=0)
        fitdat = cfp.fit(PL_log, mach_list, np.log(SNR_crit), yerr=yerr_mean/SNR_crit)
        a, b = np.exp(fitdat.popt[0]), fitdat.popt[1]
        xfit = cfp.logspace(mach_list.min(), mach_list.max(), 200)
        yfit = a * xfit**b
        label = (r'$\mathrm{SNR_{crit}}=('+cfp.round(a,2,True)+r'\pm'+cfp.round(fitdat.pstd[0],1,True)+r')\,\mathcal{M}^{'
                 +cfp.round(b,2,True)+r'}$')
        cfp.plot(x=xfit, y=yfit, ax=ax, linewidth=1.0, color='black', label=label,
                    legend_formatter=cfp.legend_formatter().add(pos=(0.25,y_of_highest_equality-i_flow_type*0.09)))
        cfp.plot(xlog=True, ylog=True, ax=ax,
                    xlabel=(r'$\mathcal{M}$' if icol==0 else None),
                    ylabel=(r'$\mathrm{SNR_{crit}}$' if icol==0 else None),
                    marker=flow_markers[flow_type],
                    label=flow_labels[flow_type])
        if icol == 0:
            ax.set_ylabel(r'$\mathrm{SNR_{crit}}$')
        #ax.set_title(flow_labels[flow_type])

        # --- row 1: PL_exp vs Mach ---
        ax = axes[1, icol]
        cfp.plot(x=mach_list, y=PL_exp_list, yerr=PL_exp_err.T, ax=ax,
                    marker=flow_markers[flow_type], linestyle='none', color=flow_colors[flow_type],
                    label=flow_labels[flow_type],
                    legend_formatter=cfp.legend_formatter().add(pos=(0.02,0.90-i_flow_type*0.09))
                    )

        # log fit code
        yerr_mean = np.abs(PL_exp_err.T).mean(axis=0)
        fitdat = cfp.fit(PL_log, mach_list, np.log(PL_exp_list), yerr=yerr_mean/PL_exp_list)
        a, b = np.exp(fitdat.popt[0]), fitdat.popt[1]
        xfit = cfp.logspace(mach_list.min(), mach_list.max(), 200)
        yfit = a * xfit**b
        label = (r'$\mathrm{\alpha}=('+cfp.round(a,2,True)+r'\pm'+cfp.round(fitdat.pstd[0],1,True)+r')\,\mathcal{M}^{'
                 +cfp.round(b,2,True)+r'}$')
        cfp.plot(x=xfit, y=yfit, ax=ax, linewidth=1.0, color='black', label=label,
                    legend_formatter=cfp.legend_formatter().add(pos=(0.25,y_of_highest_equality-i_flow_type*0.09)))
        cfp.plot(xlog=True, ylog=True, ax=ax, xlabel=r'$\mathcal{M}$')
        if icol == 0:
            ax.set_ylabel(r'$\mathrm{\alpha}$')

    fig.tight_layout()
    save = figpath + 'main/fig_4_SNR_crit_PL_exp_vs_Mach.pdf'
    fig.savefig(save, bbox_inches='tight')
    plt.close(fig)
    print("Wrote '"+save+"'.")

# =========================== FIGURE 5 ==============================================
def plot_z_proj_grid_with_noise(turb_data, flow_types=['s', 'c'], machs=[0.1, 10], snrs=[0.5, 1, 5],
                                res=256, isnr=0,  timestep=50, log=True,
                                grid_size=(9, 12), symlog=False, label_info=""):
    """
    Warning: Note that any log=/=True path is now deprecated, and likely won't work.

    plots the following grid:
                                | SNRS
                                | 0.1       1      10
    ------------------------------------------------------
    solenoidal      | Mach 0.1  |
                    | Mach 10   |
    compressive     | Mach 0.1  |
                    | Mach 10   |

    each entry is the z projection heatmap
    """

    if symlog and log:
        raise ValueError("Symlog and log can't both be true")

    cfp.load_plot_style(fontsize=20)

    j = (timestep * 3) + 2

    row_keys = [(ft, mach) for mach in machs for ft in flow_types]
    isnrs = [SNR_list_uniform.index(SNR) for SNR in snrs]

    cmap = "magma" if not symlog else 'RdBu_r'
    mach_cbars = {}

    scales = {}
    for mach in machs:
        # noisy_data[0] is the actual data in the proj

        symlog_linthresh = 1
        if symlog:
            all_z_projs = [turb_data[ft][lower_mach]['data'][isnr][j].noisy_data[0] for ft in flow_types for mach in machs]
            vmin, vmax = min(all_z_projs), max(all_z_projs)
        elif log:
            z_projs = [turb_data[ft][mach]['data'][isnr][j].noisy_data[0] for ft in flow_types for isnr in isnrs]

            flattened = np.concatenate([m.ravel() for m in z_projs])
            pos = flattened[flattened > 0]  # Note this sets all negative values to the lowest on the cbar
            #vmin, vmax = min(pos), max(pos)
            vmin, vmax = np.percentile(pos, [2, 98])
            print(f"vmin: {vmin}")
            print(f"vmax: {vmax}")
        else:
            z_projs = [turb_data[ft][mach]['data'][isnr][j].noisy_data[0] for ft in flow_types for isnr in isnrs]

            flattened = np.concatenate([m.ravel() for m in z_projs])
            #vmin, vmax = min(flattened), max(flattened)
            vmin, vmax = np.percentile(flattened, [2, 98])
            print(f"vmin: {vmin}")
            print(f"vmax: {vmax}")

        scales[mach] = (vmin, vmax)
        print(f"Mach {mach} vmin: {vmin}, vmax: {vmax}")

    fig, axes = plt.subplots(len(row_keys), len(snrs), figsize=grid_size, squeeze=False)

    for irow, (ft, mach) in enumerate(row_keys):
        vmin, vmax = scales[mach]
        for icol, isnr in enumerate(isnrs):
            ax = axes[irow, icol]
            SNR = SNR_list_uniform[isnr]

            dat = turb_data[ft][mach]['data'][isnr][j]
            if dat.SNR[0] != SNR:
                print("SNR DOESN'T MATCH")

            plot_data = dat.noisy_data[0]

            po = cfp.plot_map(plot_data, ax=ax, log=log, symlog=symlog, symlog_linthresh=symlog_linthresh,
                                axes_format=["", ""], axes_ticks=[[],[]],
                                vmin=vmin, vmax=vmax, cmap=cmap, colorbar=False)
            if mach not in mach_cbars:
                mach_cbars[mach] = ax.images[0]

            if irow == 0:
                ax.set_title(r"$\mathrm{SNR} = "+str(SNR)+r"$")
            if icol == 0:
                ax.set_ylabel(flow_labels[ft]+"\n"+r"$\mathcal{M} = "+str(mach)+r"$")

    fig.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0, hspace=0)

    save = figpath+f'main/fig_5_z_proj_grid_noise_N{res}_{label_info}.pdf'
    fig.savefig(save, bbox_inches='tight')
    plt.close(fig)
    print("Wrote '"+save+"'.", color='magenta')

    lower_vmin, lower_vmax = scales[machs[0]]
    upper_vmin, upper_vmax = scales[machs[1]]

    cfp.plot_colorbar(vmin=upper_vmin, vmax=upper_vmax, log=log, symlog=symlog, cmap=cmap, symlog_linthresh=symlog_linthresh,
                        label=r'Column density $(\mathrm{g\,cm^{-2}})$'+"\n"+r"$\mathcal{M} = "+str(machs[0])+r"$",
                        panels=3, orientation="vertical", save=figpath+f'main/fig_5_cbar_above_zero_machs.pdf')
    cfp.plot_colorbar(vmin=lower_vmin, vmax=lower_vmax,log=log,  symlog=symlog, cmap=cmap, symlog_linthresh=symlog_linthresh,
                        label=r'Column density $(\mathrm{g\,cm^{-2}})$'+"\n"+r"$\mathcal{M} = "+str(machs[1])+r"$",
                        panels=3, orientation="vertical", save=figpath+f'main/fig_5_cbar_below_zero_machs.pdf')

    load_default_plotstyle()


# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    print(datetime.now())

    # Parse the -ires argument (if it exists).
    # This is used for Parallelising, by assigning each res to a different core.
    parser = argparse.ArgumentParser()
    parser.add_argument("-ires", type=int, default=None, help="index of single res to process (for job arrays)")
    args = parser.parse_args()

    load_default_plotstyle()

    # Reduce the list of simulations to do to just the 1 per core if ires is used
    if args.ires is not None:
        resolutions = [resolutions[args.ires]]

    # reduce sim cases to load for efficency during test
    testing = False
    if testing:
        # figure 4
        flow_types = ["s", "c"]

        # Figure 5
        mach_nums = [0.1, 0.5, 2, 10]
        flow_types = ["s", "c"]

    sims = list(product(resolutions, flow_types, mach_nums))

    # load data
    turb_data = {res: {ft: {} for ft in flow_types} for res in resolutions}
    print("------------------ LOADING DATA -----------------------")
    for res, flow_type, mach in sims:
        sim_name = make_model_name_params(res, flow_type, mach)
        dillfile = f"{project_path}processed_data/{sim_name}.dill"
        with open(dillfile, "rb") as f:
            turb_data[res][flow_type][mach] = dill.load(f)
        print("READ '"+dillfile+"'.", color='green')


    # -------- FIGURES NOT TO BE USED IN PAPER --------

    """
    # Power spectrum plotted per sim
    for res, flow, mach in sims:
       plot_spectra(turb_data[flow][mach], flow, mach, res=res)


    # Machs plotted on k_noise vs SNR, plot for each flow
    for _, flow, _ in sims:
        plot_k_noise_vs_SNR_for_flows(turb_data, flow_type=flow)
    """


    # --------- FIGURES TO BE USED IN PAPER ------------

    for res in resolutions:
        print(f"------------------- RESOLUTION {res} ---------------------")

        # z projections with noise added
        print("----------------- FIGURE 5 -----------------------")
        plot_z_proj_grid_with_noise(turb_data[res],
                                    flow_types=['s', 'c'],
                                    machs=[0.5, 5],
                                    snrs=[0.3, 1, 3],
                                    symlog=False,
                                    log=False)

        """
        # Show the params of figure 3 varying.
        print("----------------- FIGURE 4-------------------------")
        plot_fit_params_grid(turb_data[res],
                            flow_types=['s', 'm', 'c'])






        # Grid of z-projections
        print("----------------- FIGURE 1 -------------------------")
        plot_z_proj_grid(turb_data[res],
                            res=res,
                            flow_types=flow_types,
                            machs=[0.1, 0.5, 2, 10],
                            grid_size=(12, 16),
                            label_info="",
                            log=True)

        # Grid of power spectra
        print("----------------- FIGURE 2 -------------------------")
        plot_spectra_grid(turb_data[res],
                            res=res,
                            flow_types=flow_types,
                            machs=[0.1, 0.5, 2, 10],
                            grid_size=(14, 16))



        # Grid of plots of x=SNR, y=knoise/kmax, across mach numbers
        print("----------------- FIGURE 3 -------------------------")
        plot_flows_grid(turb_data[res],
                            res=res,
                            flow_types=['s', 'c'],
                            grid_size=(9, 4))

        """







    print("End.")

