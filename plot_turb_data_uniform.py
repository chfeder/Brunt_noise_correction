#!/usr/bin/env python

from globals_turb_data import *
from itertools import product
from datetime import datetime
import math
import argparse
import gc



# ===================================================================

UNIVERSAL_ASYM_ERR = False # if True, all errors in graphs will be presented asymmetrically, overwriting any
                          # other option
MACH_STR = r'\mathcal{M}'

# ------------------------ HELPER FUNCTIONS --------------------------
def load_default_plotstyle():
    import matplotlib as mpl
    cfp.unload_plot_style()
    cfp.load_plot_style(fontsize=14)
    mpl.rcParams.update({
        'errorbar.capsize': 3,
        'lines.markersize': 5,
    })


def val_with_err_str(fitdat, i: int, sym=False, prec=2) -> str:
    """
    Generates strings of the form "value +/- error" for symmetrics labels, or "value^{+...}_{-...}" for
    asymetric errors.

    params:
    i: the position of the value variable in fitdat
    sym: symmetric labeling
    """

    if UNIVERSAL_ASYM_ERR:
        sym = False

    # extract variables from fitdat
    val = fitdat.popt[i]
    pstd = fitdat.pstd[i]
    perr = fitdat.perr[i]

    if val == 0:
        d = prec - 1
    else:
        d = max(prec - 1 - int(np.floor(np.log10(np.abs(val)))), 0)

    fmt = lambda x: f"{x:.{d}f}" # precision formatting

    if sym:
        return fmt(val) + r'\pm' + fmt(pstd)
    else:
        lo, hi = np.abs(perr[0]), np.abs(perr[1])
        return fmt(val) + r'^{+' + fmt(hi) + r'}_{-' + fmt(lo) + r'}'


# Generates strings such as "value +/- error" for symmetrics labels, or "value^{+...}_{-...}" for
# asymetric errors
# def val_with_err_str(val, perr, pstd, sym=True, prec=2) -> str:

#     if UNIVERSAL_ASYM_ERR:
#         sym = False

#     # the error term has the same place value as the main term
#     if val == 0:
#         d = prec - 1
#     else:
#         d = max(prec - 1 - int(np.floor(np.log10(np.abs(val)))), 0)

#     fmt = lambda x: f"{x:.{d}f}" # precision formatting

#     if sym:
#         return fmt(val) + r'\pm' + fmt(pstd)
#     else:
#         lo, hi = np.abs(perr[0]), np.abs(perr[1])
#         return fmt(val) + r'^{+' + fmt(hi) + r'}_{-' + fmt(lo) + r'}'



# ============================= FIGURE A ========================================
def plot_z_proj_grid(turb_data, res, isnr=0, flow_types=['s','m','c'], machs=[0.1, 0.5, 2, 10], timestep=50,
                        grid_size=(9, 12), log_scales={"upper_machs": True, "lower_machs": False}, label_info=""):
    # j here is the number which refers to the ret object, since they follow a pattern of
    # (0,x) -> (0,y) -> (0,z) -> (1,x) -> ...
    # for (timestep, projection LOS).

    cfp.load_plot_style(fontsize=20)

    midpt = math.floor(len(machs)/2)

    lower_machs = machs[0:midpt]
    upper_machs = machs[midpt:len(machs)]

    # set the real min and max as log boundaries
    j = (timestep * 3) + 2

    if machs == [0.1, 0.5, 2, 10]:
        dist=1.135
        lower_vmin, lower_vmax = 1/dist, dist
        dist = 60
        upper_vmin, upper_vmax = 1/dist, dist
    else:
        dist = 5 #4
        lower_vmin, lower_vmax = 1/dist, dist
        dist = 70 #65
        upper_vmin, upper_vmax = 1/dist, dist



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
                po_lower = cfp.plot_map(dat.noisy_data[0], ax=ax, log=log_scales["lower_machs"], axes_format=["", ""], axes_ticks=[[],[]],
                                    vmin=lower_vmin, vmax=lower_vmax, cmap=cmap, colorbar=False)
                ax.set_box_aspect(1)
            else:
                po_upper = cfp.plot_map(dat.noisy_data[0], ax=ax, log=log_scales["upper_machs"], axes_format=["", ""], axes_ticks=[[],[]],
                                    vmin=upper_vmin, vmax=upper_vmax, cmap=cmap, colorbar=False)
                ax.set_box_aspect(1)

            ax.text(0.5, 0.88, r"$\mathcal{M} = "+str(mach)+r"$"+"\n"+r"$\mathrm{"+str(flow_labels[ft])+r"}$",
                    ha='center', va='center', transform=ax.transAxes,
                    bbox=dict(boxstyle="round", fc="white", alpha=0.9))
            #if irow == 0:

                #ax.set_title(flow_labels[ft])

    #important: this solves the alignment problem, eliminating all the gaps
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0, hspace=0)

    if min(machs) < 1:
        save = figpath+f'main/all_machs/fig_A_z_proj_grid_'+cfp.eform(SNR,prec=1)+f'_N{res}_{label_info}.pdf'
    else:
        save = figpath+f'main/supersonic_machs/fig_A_z_proj_grid_'+cfp.eform(SNR,prec=1)+f'_N{res}_{label_info}.pdf'

    fig.savefig(save)
    plt.close(fig)
    print("Wrote '"+save+"'.", color='magenta')


    #plot colorbar
    panel_w = grid_size[0] / len(flow_types)
    panel_h = grid_size[1] / len(machs)

    if min(machs) < 1:
        higher_save = figpath+f'main/all_machs/fig_A_cbar_higher_machs.pdf'
        lower_save = figpath+f'main/all_machs/fig_A_cbar_lower_machs.pdf'
    else:
        higher_save = figpath+f'main/supersonic_machs/fig_A_cbar_higher_machs.pdf'
        lower_save = figpath+f'main/supersonic_machs/fig_A_cbar_lower_machs.pdf'

    plt.figure(figsize=(panel_w, panel_h))
    cfp.plot_colorbar(vmin=upper_vmin, vmax=upper_vmax, log=log_scales["upper_machs"], cmap=cmap, label=r'Normalised column density ($N / N_0$)',
                       panels=len(upper_machs), orientation="vertical", save=higher_save)
    plt.figure(figsize=(panel_w, panel_h))
    cfp.plot_colorbar(vmin=lower_vmin, vmax=lower_vmax, log=log_scales["lower_machs"], cmap=cmap, label=r'Normalised column density ($N / N_0$)',
                       panels=len(lower_machs), orientation="vertical", save=lower_save)


    load_default_plotstyle()


# ============================= FIGURE B ====================================
# Creates a flow_types x machs grid of power spectra plots
def plot_spectra_grid(turb_data, res, flow_types=['s', 'm', 'c'], machs=[0.1, 0.5, 2, 10], grid_size=(12, 16), skip_SNR_above_hundred=True):
    fig, axes = plt.subplots(len(machs), len(flow_types), figsize=grid_size, sharex=True, sharey=True, squeeze=False)
    for icol, ft in enumerate(flow_types):
        for irow, mach in enumerate(machs):
            plot_spectra(turb_data[ft][mach], ft, mach, ax=axes[irow, icol], save_panel=False, show_xlabel=(irow == 3),
                            show_ylabel=(icol == 0), skip_SNR_above_hundred=skip_SNR_above_hundred)

            if icol == len(flow_types) - 1:
                ax_r = axes[irow, icol].twinx()
                ax_r.set_ylabel(r"$\mathcal{M} = " + str(mach) + "$")
                ax_r.set_yticks([])

            # Add column and row titles
            if irow == 0:
                axes[irow, icol].set_title(flow_labels[ft] + rf" ($\zeta={flow_to_num[ft]}$)")
    fig.tight_layout()
    space = 0.05
    fig.subplots_adjust(wspace=space, hspace=space)

    if skip_SNR_above_hundred:
        if min(machs) < 1:
            save = figpath + f'main/all_machs/fig_B_PS_grid_N{res}.pdf'
        else:
            save = figpath + f'main/supersonic_machs/fig_B_PS_grid_N{res}.pdf'
        fig.savefig(save, bbox_inches='tight')
    else:
        fig.savefig(figpath + f'extra/PS_grid_N{res}.pdf', bbox_inches='tight')

    plt.close(fig)

def plot_spectra(turb_data, flow_type: str, mach, noise_exp=1.0, res=256,
                res_resample=1.0, normalize=False, ax=None, save_panel=True, show_xlabel=True, show_ylabel=True, skip_SNR_above_hundred=True):
    xlim = [1e-2, 2]
    ylim = [1e-13, 1e4]
    if normalize:
        ylim = [1e-8, 3e4]
    kmax = turb_data['stats'][0].avg.k.max()
    x = turb_data['stats'][0].avg.k[1:-1]/kmax
    for isnr, SNR in enumerate(SNR_list_uniform):
        pos_correct = -0.2
        if SNR > 1e2 and skip_SNR_above_hundred: continue
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
        save = figpath+'extra/fig_B_PS_N'+str(res)+flow_type+make_mach_str(mach)+noise_str+'.pdf'
    cfp.plot(xlabel=r"$k/k_\mathrm{max}$" if show_xlabel else None, ylabel=(r"$P$" if show_ylabel else None), ax=ax, xlog=True, ylog=True,
             xlim=xlim, ylim=ylim, save=save)


# =============================FIGURE C =======================================
def plot_k_noise_vs_SNR_for_flows(turb_data, flow_type, machs=mach_nums, ax=None, save_panel=True,
                                  fit_params=None, show_ylabel=True):
    # ----------------------------------------------------
    def k_noise_func(SNR, SNR_crit=1.0, PL_exp=2/3):
        #PL_exp = 2.0 / 3.0
        sharpness = 4.0
        k_noise = (1 / (1 + (SNR/SNR_crit)**(-4*PL_exp)))**(1/sharpness)
        return np.log(k_noise)
    # ----------------------------------------------------
    # similar to figure 5 (k_noise vs SNR), but for fixed 50mu voids and resolution resampling factors 0.5, 1.0, 2.0
    xlim = [5e-2, 2e4]
    mach_list = []
    SNR_crit = []
    SNR_crit_err = []
    PL_exp_list = []
    PL_exp_err = []

    # TEXT PLACEMENT CONSTANTS
    x_of_machs_labels = 0.625
    if list(machs) == list(mach_nums):
        y_of_highest_label = 0.6
    else:
        y_of_highest_label = 0.6 - 0.085 * 3

    for iMach, mach in enumerate(machs):
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
                label = r"$\mathcal{M}="+str(mach)+r"$"
                lf = cfp.legend_formatter().add(pos=(x_of_machs_labels,y_of_highest_label-iMach*0.085))
            else:
                label = None
                lf = None
            cfp.plot(x=x, y=y, yerr=np.array([yerr]).T, color=mach_colors[mach_nums.index(mach)], linestyle='none',
                        marker='o', alpha=0.8, label=label, legend_formatter=lf, ax=ax)
        snr = np.array(snr)
        k_noise = np.array(k_noise)
        k_noise_err = np.array(k_noise_err)



        # free variables to fit the function too
        params = fit_params

        fitdat = cfp.fit(k_noise_func, snr, np.log(k_noise), yerr=k_noise_err/k_noise,
                            params=params, n_random_draws=100)

        SNR_crit.append(fitdat.popt[0])
        SNR_crit_err.append(fitdat.perr[0])
        xfit = cfp.logspace(xlim[0], xlim[1], 200)
        yfit = np.exp(k_noise_func(xfit, *fitdat.popt)) # best-fit line

        # I removed label=label and legend_formatter=lf here
        cfp.plot(x=xfit, y=yfit, ax=ax, linewidth=1.0, color=mach_colors[mach_nums.index(mach)], linestyle='dotted')
    mach_list = np.array(mach_list)
    SNR_crit = np.array(SNR_crit)
    SNR_crit_err = np.array(SNR_crit_err)
    # generate final plot and save
    save = None
    if save_panel:
        save = figpath+"extra/flow_type_"+str(flow_type)+'_k_noise_vs_SNR.pdf'
    ylabel = r"$k_\mathrm{noise}/k_\mathrm{max}$" if show_ylabel else ""
    cfp.plot(xlabel=r"$\mathrm{SNR}$", ylabel=ylabel, ax=ax, save=save,
                xlog=True, ylog=True, xlim=xlim, ylim=[2e-2,2e0])


def plot_flows_grid(turb_data, res, machs=mach_nums, flow_types=['s','c'], grid_size=(12, 16), fit_params=None):
    fig, axes = plt.subplots(1, len(flow_types), figsize=grid_size, sharex=True, sharey=True, squeeze=False)
    for icol, ft in enumerate(flow_types):
            axes[0, icol].text(0.05, 0.95, flow_labels[ft] + rf" ($\zeta={flow_to_num[ft]}$)",
                           transform=axes[0, icol].transAxes,
                           ha='left', va='top')
            plot_k_noise_vs_SNR_for_flows(turb_data, ft, machs=machs, ax=axes[0, icol], save_panel=False, fit_params=fit_params, show_ylabel=True) #(icol == 0))
    fig.tight_layout()
    if min(machs) < 1:
        save = figpath + f'main/all_machs/fig_C_flows_grid_N{res}.pdf'
    else:
        save = figpath + f'main/supersonic_machs/fig_C_flows_grid_N{res}.pdf'
    fig.savefig(save, bbox_inches='tight')
    plt.close(fig)



# =========================== FIGURE D ==============================================
# somehwat of a copy of figure C's, except this function doesn't plot, it only fits the data.
# used to create the subpanels.
def fit_k_noise_vs_SNR(turb_data, flow_type, machs, fit_params):
    # ----------------------------------------------------
    def k_noise_func(SNR, SNR_crit=1.0, PL_exp=2/3):
        sharpness = 4.0
        k_noise = (1 / (1 + (SNR/SNR_crit)**(-4*PL_exp)))**(1/sharpness)
        return np.log(k_noise)
    # ----------------------------------------------------
    mach_list = []
    SNR_crit = []
    SNR_crit_err = []
    PL_exp_list = []
    PL_exp_err = []

    for iMach, mach in enumerate(machs):
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
        params = fit_params

        fitdat = cfp.fit(k_noise_func, snr, np.log(k_noise), yerr=k_noise_err/k_noise,
                         params=params, n_random_draws=100)

        SNR_crit.append(fitdat.popt[0]); SNR_crit_err.append(fitdat.perr[0])
        PL_exp_list.append(fitdat.popt[1]); PL_exp_err.append(fitdat.perr[1])

    return (np.array(mach_list), np.array(SNR_crit), np.array(SNR_crit_err),
            np.array(PL_exp_list), np.array(PL_exp_err))

def plot_fit_params_grid(turb_data, flow_types=['s', 'm', 'c'], machs=mach_nums, fit_params=None):
    """
    2-panel graph as follows:


    SNR_crit        |   1   |
    PL_exp = lambda |   2   |
    ---------------------------------
                       Mach

    each panel has the param as the y-axis, and all have the mach as x-axis

    """
    fig, axes = plt.subplots(2, 1, figsize=(6, 8), squeeze=False, sharex=True)


    #def PL(x, a, b):
    #    return a*x**b

    def PL_log(x, a, b):
        return np.log(a) + b * np.log(x)


    fit_params_ab_SNR_crit = {
        'a': [1e-3, 100, 1e3],
        'b': [-10, -1, 10],
    }

    fit_params_ab_alpha = {
        'a': [1e-3, 1, 1e3],
        'b': [-10, 0.1, 10],
    }

    # TEXT PLACEMENT CONSTANTS
    y_of_top_label = 0.9

    x_of_SNR_crit_eq_label = 0.0
    x_of_alpha_eq_label = 0


    if np.array_equal(machs, mach_nums):      # full set of Machs
        ylim_SNR_crit = [4, 8e4]
        ylim_alpha    = [0.25, 1.5]
    else:
        ylim_SNR_crit = [4, 5e3]
        ylim_alpha    = [0.4, 1.1]


    for i_flow_type, flow_type in enumerate(flow_types):
        mach_list, SNR_crit, SNR_crit_err, PL_exp_list, PL_exp_err = fit_k_noise_vs_SNR(turb_data, flow_type, machs, fit_params)

        mach_list, mach_std = get_true_mach(flow_type, mach_list) #update with true mach values

        icol = 0

        # --- row 0: SNR_crit vs Mach ---
        ax = axes[0, icol]
        cfp.plot(x=mach_list, y=SNR_crit, xerr=mach_std, yerr=SNR_crit_err.T, ax=ax,
                    marker='o', linestyle=':', color=flow_colors[flow_type])

        yerr_mean = np.abs(SNR_crit_err.T).mean(axis=0)
        fitdat = cfp.fit(PL_log, mach_list, np.log(SNR_crit), xerr=mach_std, yerr=yerr_mean/SNR_crit, params=fit_params_ab_SNR_crit)
        a, b = fitdat.popt[0], fitdat.popt[1]
        xfit = cfp.logspace(mach_list.min(), mach_list.max(), 200)
        yfit = a * xfit**b

        zeta_str = f'{flow_to_num[flow_type]:.1f}'
        label = (r'$\zeta = ' + zeta_str + r';\ \mathrm{SNR_{crit}}=(' + val_with_err_str(fitdat, 0, sym=True) + r')\,'
                    + MACH_STR + r'^{(' + val_with_err_str(fitdat, 1, sym=True) + r')}$')

        cfp.plot(x=xfit, y=yfit, ax=ax, linewidth=1.0, label=label,
                    legend_formatter=cfp.legend_formatter().add(pos=(x_of_SNR_crit_eq_label,y_of_top_label-i_flow_type*0.09)),
                    color=flow_colors[flow_type], linestyle='-',
                    ylim=ylim_SNR_crit,
                    )
        cfp.plot(xlog=True, ylog=True, ax=ax,
                    xlabel=None,
                    ylabel=(r'$\mathrm{SNR_{crit}}$' if icol==0 else None),
                    marker='o',)
        if icol == 0:
            ax.set_ylabel(r'$\mathrm{SNR_{crit}}$')
        #ax.set_title(flow_labels[flow_type])

        # --- row 1: PL_exp vs Mach ---
        ax = axes[1, icol]
        cfp.plot(x=mach_list, y=PL_exp_list, xerr=mach_std, yerr=PL_exp_err.T, ax=ax,
                    marker='o', linestyle=':', color=flow_colors[flow_type])


        yerr_mean = np.abs(PL_exp_err.T).mean(axis=0)
        fitdat = cfp.fit(PL_log, mach_list, np.log(PL_exp_list), xerr=mach_std, yerr=yerr_mean/PL_exp_list, params=fit_params_ab_alpha)
        a, b = fitdat.popt[0], fitdat.popt[1]
        xfit = cfp.logspace(mach_list.min(), mach_list.max(), 200)
        yfit = a * xfit**b
        label = (r'$\zeta = ' + zeta_str + r';\ \mathrm{\lambda}=(' + val_with_err_str(fitdat, 0, sym=True) + r')\,'
                    + MACH_STR + r'^{(' + val_with_err_str(fitdat, 1, sym=True) + r')}$')
        cfp.plot(x=xfit, y=yfit, ax=ax, linewidth=1.0,label=label,
                    legend_formatter=cfp.legend_formatter().add(pos=(x_of_alpha_eq_label,y_of_top_label-i_flow_type*0.09)),
                    color=flow_colors[flow_type], linestyle='-',
                    ylim=ylim_alpha,)
        cfp.plot(xlog=True, ylog=True, ax=ax, xlabel=r'$\mathcal{M}$')
        if icol == 0:
            ax.set_ylabel(r'$\mathrm{\lambda}$')

    axes[0, 0].set_xlabel('')
    axes[0, 0].tick_params(axis='x', labelbottom=False)
    fig.tight_layout()
    if min(machs) < 1:
        save = figpath + 'main/all_machs/fig_D_SNR_crit_PL_exp_vs_Mach.pdf'
    else:
        save = figpath + 'main/supersonic_machs/fig_D_SNR_crit_PL_exp_vs_Mach.pdf'
    fig.savefig(save, bbox_inches='tight')
    plt.close(fig)
    print("Wrote '"+save+"'.")

# =========================== FIGURE E ==============================================
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

    scales = {
        0.2: (0.955, 1.045), #tailoured
        1: (0, 2),
        5.0: (-3, 5), #tailoured
        10: (-6, 8),
    }

    symlog_linthresh = 1


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

    if min(machs) < 1:
        save = figpath+f'main/all_machs/fig_E_z_proj_grid_noise_N{res}_{label_info}.pdf'
    else:
        save = figpath+f'main/supersonic_machs/fig_E_z_proj_grid_noise_N{res}_{label_info}.pdf'
    fig.savefig(save)
    plt.close(fig)
    print("Wrote '"+save+"'.", color='magenta')

    lower_vmin, lower_vmax = scales[machs[0]]
    upper_vmin, upper_vmax = scales[machs[1]]

    panel_w = grid_size[0] / len(snrs)
    panel_h = grid_size[1] / len(row_keys)

    if min(machs) < 1:
        save_upper = figpath+f'main/all_machs/fig_E_cbar_higher_machs.pdf'
        save_lower = figpath+f'main/all_machs/fig_E_cbar_lower_machs.pdf'
    else:
        save_upper = figpath+f'main/supersonic_machs/fig_E_cbar_higher_machs.pdf'
        save_lower = figpath+f'main/supersonic_machs/fig_E_cbar_lower_machs.pdf'

    plt.figure(figsize=(panel_w, panel_h))
    cfp.plot_colorbar(vmin=upper_vmin, vmax=upper_vmax, log=log, symlog=symlog, cmap=cmap, symlog_linthresh=symlog_linthresh,
                        label=r'Normalised column density ($N / N_0$)',
                        panels=2, orientation="vertical", save=save_upper)
    plt.figure(figsize=(panel_w, panel_h))
    cfp.plot_colorbar(vmin=lower_vmin, vmax=lower_vmax,log=log,  symlog=symlog, cmap=cmap, symlog_linthresh=symlog_linthresh,
                        label=r'Normalised column density ($N / N_0$)',
                        panels=2, orientation="vertical", save=save_lower)

    load_default_plotstyle()




# =============================FIGURE F =======================================
"""
Mach 0.1 vs Mach 10

"""


def plot_k_noise_vs_SNR_for_machs(turb_data, mach, ax=None, save_panel=True, fit_params=None):
    # ----------------------------------------------------
    def k_noise_func(SNR, SNR_crit=1.0, PL_exp=2/3):
        #PL_exp = 2.0 / 3.0
        sharpness = 4.0
        k_noise = (1 / (1 + (SNR/SNR_crit)**(-4*PL_exp)))**(1/sharpness)
        return np.log(k_noise)
    # ----------------------------------------------------

    # TEXT PLACEMENT CONSTANTS
    x_of_mach_labels = 0.625 #formerly 0.65
    y_of_heighest_mach_label = 0.7#formerly 0.65


    # similar to figure 5 (k_noise vs SNR), but for fixed 50mu voids and resolution resampling factors 0.5, 1.0, 2.0
    xlim = [5e-2, 2e4]
    zeta_list = []
    SNR_crit = []
    SNR_crit_err = []
    PL_exp_list = []
    PL_exp_err = []

    for iZeta, flow_type in enumerate(flow_types):
        kmax = turb_data[flow_type][mach]['stats'][0].avg.k.max() #this should be context for the uniform experiemtn

        print(f"kmax is {kmax}")

        zeta_list.append(flow_type)
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
                label = r"$\mathcal{\zeta}="+str(flow_to_num[flow_type])+r"$"
                lf = cfp.legend_formatter().add(pos=(x_of_mach_labels,y_of_heighest_mach_label-iZeta*0.085))
            else:
                label = None
                lf = None
            cfp.plot(x=x, y=y, yerr=np.array([yerr]).T, color=zeta_colours[iZeta], linestyle='none',
                        marker='o', alpha=0.8, label=label, legend_formatter=lf, ax=ax)
        snr = np.array(snr)
        k_noise = np.array(k_noise)
        k_noise_err = np.array(k_noise_err)

        # free variables to fit the function too
        # format: [lower bound, starting guess, upper bound]
        params = fit_params
        fitdat = cfp.fit(k_noise_func, snr, np.log(k_noise), yerr=k_noise_err/k_noise,
                            params=params, n_random_draws=100)

        SNR_crit.append(fitdat.popt[0])
        SNR_crit_err.append(fitdat.perr[0])
        xfit = cfp.logspace(xlim[0], xlim[1], 200)
        yfit = np.exp(k_noise_func(xfit, *fitdat.popt)) # best-fit line

        # I removed label=label and legend_formatter=lf here
        cfp.plot(x=xfit, y=yfit, ax=ax, linewidth=1.0, color=zeta_colours[iZeta], linestyle='dotted')
    zeta_list = np.array(zeta_list)
    SNR_crit = np.array(SNR_crit)
    SNR_crit_err = np.array(SNR_crit_err)
    # generate final plot and save
    save = None
    if save_panel:
        save = figpath+"extra/flow_type_"+str(flow_type)+'_k_noise_vs_SNR.pdf'
    cfp.plot(xlabel=r"$\mathrm{SNR}$", ylabel=r"$k_\mathrm{noise}/k_\mathrm{max}$", ax=ax, save=save,
                xlog=True, ylog=True, xlim=xlim, ylim=[2e-2,2e0])


def plot_machs_grid(turb_data, res, machs=[0.1, 10], grid_size=(12, 16), fit_params=None):
    fig, axes = plt.subplots(1, len(machs), figsize=grid_size, sharex=True, sharey=True, squeeze=False)
    for icol, mach in enumerate(machs):
            axes[0, icol].text(0.05, 0.95, r'$\mathcal{M}=' + str(mach) + r'$',
                           transform=axes[0, icol].transAxes,
                           ha='left', va='top')
            plot_k_noise_vs_SNR_for_machs(turb_data, mach, ax=axes[0, icol], save_panel=False, fit_params=fit_params)
    fig.tight_layout()
    if min(machs) < 1:
        save = figpath + f'main/all_machs/fig_F_flows_grid_N{res}.pdf'
    else:
        save = figpath + f'main/supersonic_machs/fig_F_flows_grid_N{res}.pdf'
    fig.savefig(save, bbox_inches='tight')


# =========================== FIGURE G ==============================================
# somehwat of a copy of figure F's, except this function doesn't plot, it only fits the data.
# used to create the subpanels.
def fit_k_noise_vs_SNR_for_zeta(turb_data, mach, fit_params):
    # ----------------------------------------------------
    def k_noise_func(SNR, SNR_crit=1.0, PL_exp=2/3):
        sharpness = 4.0
        k_noise = (1 / (1 + (SNR/SNR_crit)**(-4*PL_exp)))**(1/sharpness)
        return np.log(k_noise)
    # ----------------------------------------------------
    flows_list = []
    SNR_crit = []
    SNR_crit_err = []
    PL_exp_list = []
    PL_exp_err = []

    for iFlow, flow_type in enumerate(flow_types):
        kmax = turb_data[flow_type][mach]['stats'][0].avg.k.max()
        flows_list.append(flow_type)
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
        params = fit_params
        fitdat = cfp.fit(k_noise_func, snr, np.log(k_noise), yerr=k_noise_err/k_noise,
                         params=params, n_random_draws=100)

        SNR_crit.append(fitdat.popt[0]); SNR_crit_err.append(fitdat.perr[0])
        PL_exp_list.append(fitdat.popt[1]); PL_exp_err.append(fitdat.perr[1])

    return (np.array(flows_list), np.array(SNR_crit), np.array(SNR_crit_err),
            np.array(PL_exp_list), np.array(PL_exp_err))

def plot_fit_params_grid_zeta(turb_data, machs=[0.1, 1, 10], fit_params=None):
    """
    2-panel graph as follows:


    SNR_crit        |   1   |
    PL_exp = lambda |   2   |
    ---------------------------------
                       Zeta

    each panel has the param as the y-axis, and all have the zets as x-axis

    the xaxis is on a linear scale, the y-axis is on a log scale.

    """
    fig, axes = plt.subplots(2, 1, figsize=(6, 8), squeeze=False, sharex=True)

    def exp_func(x, b, c):
        return np.exp(b * (1 - x)) + c

    def log_exp_func(x, b, c):
        return np.log(np.exp(b * (1 - x)) + c)

    def lin_func(x, a, b):
       return a + x * b

    # Text placement constants
    x_of_upper_row_flow_type = 0


    label_machs = [0.1, 1, 10]
    plotting_all_machs = (list(machs) == list(mach_nums))
    iLabel = 0

    if plotting_all_machs:
        ylim_SNR_crit = [5, 1e5]
        ylim_alpha    = [0.29, 1.4]
    else:
        ylim_SNR_crit = [3, 2e4]
        ylim_alpha    = [0.4, 1.5]


    for iMach, mach in enumerate(machs):
        show_label = (mach in label_machs) if plotting_all_machs else True

        flows_list, SNR_crit, SNR_crit_err, PL_exp_list, PL_exp_err = fit_k_noise_vs_SNR_for_zeta(turb_data, mach, fit_params=fit_params)

        icol = 0

        # --- row 0: SNR_crit vs Mach ---
        ax = axes[0, icol]

        ylim = ylim_SNR_crit

        yerr_mean = np.abs(SNR_crit_err.T).mean(axis=0)

        # format: [lower bound, starting guess, upper bound]
        exp_fit_params = {
            'b': [0, 5, 10],
            'c': [0, 10, 1e4],
        }
        fitdat = cfp.fit(log_exp_func, flow_types_nums, np.log(SNR_crit), yerr=yerr_mean/SNR_crit, params=exp_fit_params)
        b, c = fitdat.popt[0], fitdat.popt[1]
        xfit = np.linspace(0, 1, 200)
        yfit = exp_func(xfit, b, c)
        label = (r'$\mathcal{M} = ' + str(mach) + r';\; \mathrm{SNR_{crit}} = (' + val_with_err_str(fitdat, 1, sym=True)
         + r') + e^{(' + val_with_err_str(fitdat, 0, sym=True)
         + r')(1-\zeta)}$')

        if show_label:
            cfp.plot(x=flow_types_nums, y=SNR_crit, yerr=SNR_crit_err.T, ax=ax,
                        marker='o', linestyle='dotted', color=mach_colors[mach_nums.index(mach)],
                        label=label,
                        legend_formatter=cfp.legend_formatter().add(pos=(0.02,0.90-iLabel*0.09)),
                        ylim=ylim
                        )
            cfp.plot(x=xfit, y=yfit, ax=ax, linewidth=1.0, color=mach_colors[mach_nums.index(mach)], label=None)
            cfp.plot(xlog=False, ylog=True, ax=ax,
                        xlabel=None,
                        ylabel=(r'$\mathrm{SNR_{crit}}$' if icol==0 else None),
                        marker='o')
        else:
            cfp.plot(x=flow_types_nums, y=SNR_crit, yerr=SNR_crit_err.T, ax=ax,
                        marker='o', linestyle='dotted', color=mach_colors[mach_nums.index(mach)],
                        ylim=ylim
                        )
            cfp.plot(x=xfit, y=yfit, ax=ax, linewidth=1.0, color=mach_colors[mach_nums.index(mach)], label=None)
            cfp.plot(xlog=False, ylog=True, ax=ax,
                        xlabel=None,
                        ylabel=(r'$\mathrm{SNR_{crit}}$' if icol==0 else None),
                        marker='o',)
        if icol == 0:
            ax.set_ylabel(r'$\mathrm{SNR_{crit}}$')
        #ax.set_title(flow_labels[flow_type])

        # --- row 1: PL_exp vs Mach ---

        ax = axes[1, icol]

        ylim = ylim_alpha

        yerr_mean = np.abs(PL_exp_err.T).mean(axis=0)
        fitdat = cfp.fit(lin_func, flow_types_nums, PL_exp_list, yerr=yerr_mean)
        a, b = fitdat.popt[0], fitdat.popt[1]
        xfit = np.linspace(0, 1, 200)
        yfit = lin_func(xfit, a, b)



        label = (r'$\mathcal{M} = ' + str(mach) + r';\; \lambda = (' + val_with_err_str(fitdat, 0, sym=True)
                + r') + (' + val_with_err_str(fitdat, 1, sym=True) + r')\,\zeta$')
        if show_label:
            cfp.plot(x=flow_types_nums, y=PL_exp_list, yerr=PL_exp_err.T, ax=ax,
                        marker='o', linestyle='dotted', color=mach_colors[mach_nums.index(mach)],
                        label=label,
                        legend_formatter=cfp.legend_formatter().add(pos=(0.02,0.90-iLabel*0.09)),
                        ylim=ylim
                        )
            cfp.plot(x=xfit, y=yfit, ax=ax, linewidth=1.0, color=mach_colors[mach_nums.index(mach)], label=None)
        else:
            cfp.plot(x=flow_types_nums, y=PL_exp_list, yerr=PL_exp_err.T, ax=ax,
                        marker='o', linestyle='dotted', color=mach_colors[mach_nums.index(mach)],
                        ylim=ylim
                        )
            cfp.plot(x=xfit, y=yfit, ax=ax, linewidth=1.0, color=mach_colors[mach_nums.index(mach)])
        cfp.plot(xlog=False, ylog=True, ax=ax, xlabel=r'$\zeta$')
        if icol == 0:
            ax.set_ylabel(r'$\mathrm{\lambda}$')

        if show_label:
            iLabel += 1

    axes[0, 0].set_xlabel('')
    axes[0, 0].tick_params(axis='x', labelbottom=False)
    fig.tight_layout()
    if min(machs) < 1:
        save = figpath + 'main/all_machs/fig_G_SNR_crit_and_alpha_vs_zeta.pdf'
    else:
        save = figpath + 'main/supersonic_machs/fig_G_SNR_crit_and_alpha_vs_zeta.pdf'


    fig.savefig(save, bbox_inches='tight')
    plt.close(fig)
    print("Wrote '"+save+"'.")


# =========================== FIGURE H ==============================================
def plot_fig_H(turb_data, machs=[0.1, 1, 10], fit_params=None):
    """
    2-panel graph:

    a       |   1   |

    b       |   2   |

    -----------------
              Mach
    """

    x_of_lam_eq = 0.25
    x_of_eta_eq = x_of_lam_eq + 0.2

    fig, axes = plt.subplots(2, 1, figsize=(6, 8), squeeze=False, sharex=True)

    def log_exp_func(x, b, c):
        return np.log(np.exp(b * (1 - x)) + c)

    def log_func(x, a, b):
        return a * np.log(x) + b

    def PL_log(x, A, p):
        return np.log(A) + p * np.log(x)

    lam_list, lam_perr_list, eta_list, eta_perr_list = [], [], [], []
    true_mach_list, true_mach_std_list = [], []


    # For each mach, fit SNR_crit(zeta) then plot
    for iMach, mach in enumerate(machs):
        flows_list, SNR_crit, SNR_crit_err, PL_exp_list, PL_exp_err = fit_k_noise_vs_SNR_for_zeta(turb_data, mach, fit_params=fit_params)

        true_mach, true_mach_std = get_true_mach_agg_over_flows(flows_list, mach)
        true_mach_list.append(true_mach)
        true_mach_std_list.append(true_mach_std)

        yerr_mean = np.abs(SNR_crit_err.T).mean(axis=0)

        exp_fit_params = {
            'b': [0, 5, 10],
            'c': [0, 10, 1e4],
        }
        fitdat = cfp.fit(log_exp_func, flow_types_nums, np.log(SNR_crit), yerr=yerr_mean/SNR_crit, params=exp_fit_params)

        lam_list.append(fitdat.popt[0])
        lam_perr_list.append(fitdat.perr[0])
        eta_list.append(fitdat.popt[1])
        eta_perr_list.append(fitdat.perr[1])

        color = mach_colors[mach_nums.index(mach)]
        cfp.plot(x=true_mach, y=fitdat.popt[0], xerr=true_mach_std, yerr=np.array([np.abs(fitdat.perr[0])]).T, ax=axes[1, 0],
                 color=color, marker='o', linestyle='none')
        cfp.plot(x=true_mach, y=fitdat.popt[1], xerr=true_mach_std, yerr=np.array([np.abs(fitdat.perr[1])]).T, ax=axes[0, 0],
                 color=color, marker='o', linestyle='none')

    #machs_arr = np.array(machs)
    machs_arr = np.array(true_mach_list)
    machs_std_arr = np.array(true_mach_std_list)
    lam_list = np.array(lam_list)
    lam_perr_list = np.array(lam_perr_list)
    eta_list = np.array(eta_list)
    eta_perr_list = np.array(eta_perr_list)
    order = np.argsort(machs_arr)


    xfit = cfp.logspace(min(machs), max(machs), 200)

    # FITTING LINES OF BEST FIT
    # row 1: lambda vs Mac
    ax = axes[1, 0]
    ax.plot(machs_arr[order], lam_list[order], linestyle=':', color='k', linewidth=1.0, zorder=2)

    # format: [lower bound, starting guess, upper bound]
    log_fit_params = {
        'a': [-100, -1, 0],
        'b': [-100, 1, 100],
    }
    yerr_mean = np.abs(lam_perr_list).mean(axis=1)
    fitdat = cfp.fit(log_func, machs_arr, lam_list, xerr=machs_std_arr, yerr=yerr_mean, params=log_fit_params)
    label = (r'$b = (' + val_with_err_str(fitdat, 0, sym=True)
            + r')\ln' + MACH_STR + r' + (' + val_with_err_str(fitdat, 1, sym=True) + r')$')

    cfp.plot(x=xfit, y=log_func(xfit, *fitdat.popt), ax=ax, linewidth=1.0, color='k', label=label,
             legend_formatter=cfp.legend_formatter().add(pos=(x_of_lam_eq, 0.90)))
    cfp.plot(xlog=True, ylog=False, ax=ax, xlabel=r'$\mathcal{M}$')
    ax.set_ylabel(r'$b$')

    # row 0: eta vs Mach
    ax = axes[0, 0]
    ax.plot(machs_arr[order], eta_list[order], linestyle=':', color='k', linewidth=1.0, zorder=2)

    pl_fit_params = {
        'A': [1e-6, 1, 1e5],
        'p': [-5, 0.5, 5],
    }
    yerr_mean = np.abs(eta_perr_list).mean(axis=1)
    fitdat = cfp.fit(PL_log, machs_arr, np.log(eta_list), xerr=machs_std_arr, yerr=yerr_mean/eta_list, params=pl_fit_params)
    A, p = fitdat.popt[0], fitdat.popt[1]
    label = (r'$a = (' + val_with_err_str(fitdat, 0, sym=True)
            + r')\,' + MACH_STR + r'^{(' + val_with_err_str(fitdat, 1, sym=True) + r')}$')
    cfp.plot(x=xfit, y=A * xfit**p, ax=ax, linewidth=1.0, color='k', label=label,
             legend_formatter=cfp.legend_formatter().add(pos=(x_of_eta_eq, 0.90)))
    cfp.plot(xlog=True, ylog=True, ax=ax, xlabel=None)
    ax.set_ylabel(r'$a$')

    axes[0, 0].set_xlabel('')
    axes[0, 0].tick_params(axis='x', labelbottom=False)
    fig.tight_layout()
    if min(machs) < 1:
        save = figpath + 'main/all_machs/fig_H_SNR_crit_fit_params_vs_mach.pdf'
    else:
        save = figpath + 'main/supersonic_machs/fig_H_SNR_crit_fit_params_vs_mach.pdf'

    fig.savefig(save, bbox_inches='tight')
    plt.close(fig)
    print("Wrote '"+save+"'.")


# =========================== FIGURE I ==============================================
def plot_fig_I(turb_data, machs=[0.1, 1, 10], fit_params=None):
    """
    2-panel graph:

    alpha_0  |   1   |
    alpha_1  |   2   |
    ------------------
                Mach
    """
    fig, axes = plt.subplots(2, 1, figsize=(6, 8), squeeze=False, sharex=True)

    def lin_func(x, a, b):
        return a + x * b

    def log_func(x, a, b):
        return a * np.log(x) + b

    def log_func_no_intercept(x, a):
        return a * np.log(x)

    a0_list, a0_perr_list, a1_list, a1_perr_list = [], [], [], []
    true_mach_list, true_mach_std_list = [], []

    # For each mach, fit alpha(zeta) then plot
    for iMach, mach in enumerate(machs):
        flows_list, SNR_crit, SNR_crit_err, PL_exp_list, PL_exp_err = fit_k_noise_vs_SNR_for_zeta(turb_data, mach, fit_params=fit_params)

        true_mach, true_mach_std = get_true_mach_agg_over_flows(flows_list, mach)
        true_mach_list.append(true_mach)
        true_mach_std_list.append(true_mach_std)

        yerr_mean = np.abs(PL_exp_err.T).mean(axis=0)
        fitdat = cfp.fit(lin_func, flow_types_nums, PL_exp_list, yerr=yerr_mean)

        a0_list.append(fitdat.popt[0])
        a0_perr_list.append(fitdat.perr[0])
        a1_list.append(fitdat.popt[1])
        a1_perr_list.append(fitdat.perr[1])

        color = mach_colors[mach_nums.index(mach)]
        cfp.plot(x=true_mach, y=fitdat.popt[0], xerr=true_mach_std, yerr=np.array([np.abs(fitdat.perr[0])]).T, ax=axes[0, 0],
                 color=color, marker='o', linestyle='none')
        cfp.plot(x=true_mach, y=fitdat.popt[1], xerr=true_mach_std, yerr=np.array([np.abs(fitdat.perr[1])]).T, ax=axes[1, 0],
                 color=color, marker='o', linestyle='none')

    machs_arr = np.array(true_mach_list)
    machs_std_arr = np.array(true_mach_std_list)
    a0_list = np.array(a0_list)
    a0_perr_list = np.array(a0_perr_list)
    a1_list = np.array(a1_list)
    a1_perr_list = np.array(a1_perr_list)
    order = np.argsort(machs_arr)

    xfit = cfp.logspace(min(machs), max(machs), 200)

    log_fit_params = {
        'a': [-100, 0.1, 100],
        'b': [-100, 0.1, 100],
    }
    log_fit_params_no_intercept = {
        'a': [-100, 0.1, 100],
    }

    # d keeps its intercept only when fitting the full set of Machs
    d_has_intercept = list(machs) == list(mach_nums)

    d_func = log_func if d_has_intercept else log_func_no_intercept
    d_params = log_fit_params if d_has_intercept else log_fit_params_no_intercept

    # (axis, data, errors, symbol, fit function, fit params, has intercept)
    panels = [
        (axes[0, 0], a0_list, a0_perr_list, r'c', log_func, log_fit_params, True),
        (axes[1, 0], a1_list, a1_perr_list, r'd', d_func, d_params, d_has_intercept),
    ]

    for ax, y_list, y_perr_list, sym, func, params, has_intercept in panels:
        ax.plot(machs_arr[order], y_list[order], linestyle=':', color='k', linewidth=1.0, zorder=2)

        yerr_mean = np.abs(y_perr_list).mean(axis=1)
        fitdat = cfp.fit(func, machs_arr, y_list, xerr=machs_std_arr, yerr=yerr_mean, params=params)

        a_str = val_with_err_str(fitdat, 0, sym=True)

        label = r'$' + sym + r' = (' + a_str + r')\ln' + MACH_STR
        if has_intercept:
            b_str = val_with_err_str(fitdat, 1, sym=True)
            label += r' + (' + b_str + r')'
        label += r'$'

        cfp.plot(x=xfit, y=func(xfit, *fitdat.popt), ax=ax, linewidth=1.0, color='k', label=label,
                 legend_formatter=cfp.legend_formatter().add(pos=(0.02, 0.90)))

        cfp.plot(xlog=True, ylog=False, ax=ax, xlabel=r'$\mathcal{M}$')
        ax.set_ylabel(r'$' + sym + r'$')

    axes[0, 0].set_xlabel('')
    axes[0, 0].tick_params(axis='x', labelbottom=False)
    fig.tight_layout()
    if min(machs) < 1:
        save = figpath + 'main/all_machs/fig_I_alpha_fit_params_vs_mach.pdf'
    else:
        save = figpath + 'main/supersonic_machs/fig_I_alpha_fit_params_vs_mach.pdf'

    fig.savefig(save, bbox_inches='tight')
    plt.close(fig)
    print("Wrote '"+save+"'.")


# ========================DATA LOADING HELPER FUNCTIONS ==================================


def load_data(resolutions, flow_types, machs):
    sims = list(product(resolutions, flow_types, machs))
    turb_data = {res: {ft: {} for ft in flow_types} for res in resolutions}

    print("------------------ LOADING DATA -----------------------")
    for res, flow_type, mach in sims:
        sim_name = make_model_name_params(res, flow_type, mach)
        dillfile = f"{project_path}processed_data/{sim_name}.dill"
        with open(dillfile, "rb") as f:
            turb_data[res][flow_type][mach] = dill.load(f)
        print("READ '"+dillfile+"'.", color='green')

    return turb_data # this is ok in runtime since python returns by reference


def load_data_stats_only(resolutions, flow_types, machs):
    turb_data = {res: {ft: {} for ft in flow_types} for res in resolutions}

    print("------------------ LOADING DATA (stats only) -----------------------")
    for res, flow_type, mach in product(resolutions, flow_types, machs):
        sim_name = make_model_name_params(res, flow_type, mach)
        dillfile = f"{project_path}processed_data/{sim_name}.dill"
        with open(dillfile, "rb") as f:
            full = dill.load(f)
        turb_data[res][flow_type][mach] = {'stats': full['stats']}
        del full
        gc.collect()
        print("READ '"+dillfile+"'.", color='green')

    return turb_data



# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    # ================ SETTINGS ================

    # Defines the search space for cfp.fit calls
    # format: [lower bound, starting guess, upper bound]
    fit_params = {
        'SNR_crit': [0.001, 10, 1e4],   #original: [0.1, 1.0, 50]
        'PL_exp':   [0.001, 1.0, 100]  #originally fixed to 2/3
    }


    # Determine which Mach range plots to show
    full_mach_range = True
    supersonic_mach_range = True

    # ============= END OF SETTINGS ==============

    print(datetime.now())

    # Parse the -ires argument (if it exists).
    # This is used for Parallelising, by assigning each res to a different core.
    parser = argparse.ArgumentParser()
    parser.add_argument("-ires", type=int, default=None, help="index of single res to process (for job arrays)")
    args = parser.parse_args()

    load_default_plotstyle()


    # Decide which subset of the data to load
    # Loading all data without parrellisation causes a crash, figures are divided into two subcases, loaded seperately
    if args.ires is not None:
        if args.ires == 0:
            data_and_stats_plots = True
            stats_only_plots = False
        elif args.ires == 1:
            data_and_stats_plots = False
            stats_only_plots = True
        else:
            print(f"Warning: args.ires is {args.ires}, which may lead to undefined behaviour", color=red)
    else:
        data_and_stats_plots = True
        stats_only_plots = True


    if data_and_stats_plots:
        # load data partially
        print("Loading SMC Data", color='blue')
        turb_data = load_data(
            resolutions=resolutions,
            flow_types=["s", 'm', "c"],
            machs=(mach_nums if full_mach_range else supersonic_mach_nums),
        )

        for res in resolutions:
            print(f"------------------- RESOLUTION {res} ---------------------")

            if full_mach_range:
                print(f"----------------- PLOTTING FULL MACH RANGE FIGS ---------------------",
                        color='yellow')

                # Grid of z-projections
                print("----------------- FIGURE A -------------------------")
                plot_z_proj_grid(turb_data[res],
                                    res=res,
                                    flow_types=['s', 'm', 'c'],
                                    machs=[0.1, 0.5, 2, 10],
                                    grid_size=(12, 16),
                                    label_info="",
                                    log_scales={
                                        "upper_machs": True,
                                        "lower_machs": True
                                    },)



                # z projections with noise added
                print("----------------- FIGURE E -----------------------")
                plot_z_proj_grid_with_noise(turb_data[res],
                                            flow_types=['s', 'c'],
                                            machs=[0.2, 5],
                                            snrs=[0.3, 1, 3],
                                            symlog=False,
                                            log=False)

            if supersonic_mach_range:
                print(f"----------------- PLOTTING SUPERSONIC MACH RANGE FIGS ---------------------",
                        color='yellow')

                # Grid of z-projections
                print("----------------- FIGURE A -------------------------")
                plot_z_proj_grid(turb_data[res],
                                    res=res,
                                    flow_types=['s', 'm', 'c'],
                                    machs=[1, 2, 5, 10],
                                    grid_size=(12, 16),
                                    label_info="",
                                    log_scales={
                                        "upper_machs": True,
                                        "lower_machs": True
                                    },)



                # z projections with noise added
                print("----------------- FIGURE E -----------------------")
                plot_z_proj_grid_with_noise(turb_data[res],
                                            flow_types=['s', 'c'],
                                            machs=[1, 10],
                                            snrs=[0.3, 1, 3],
                                            symlog=False,
                                            log=False)


            # deallocate data
            # relevant only if zeta is also called later
            del turb_data
            gc.collect()




    if stats_only_plots:
        print("Loading Zeta Data", color='blue')
        # load data partially

        turb_data = load_data_stats_only(
            resolutions=resolutions,
            flow_types=flow_types,
            machs=(mach_nums if full_mach_range else supersonic_mach_nums),
        )


        for res in resolutions:
            print(f"------------------- RESOLUTION {res} ---------------------")


            if full_mach_range:
                print(f"----------------- PLOTTING FULL MACH RANGE FIGS ---------------------",
                        color='yellow')

                # Grid of power spectra
                print("----------------- FIGURE B -------------------------")
                plot_spectra_grid(turb_data[res],
                                    res=res,
                                    flow_types=['s', 'm', 'c'],
                                    machs=[0.1, 0.5, 2, 10],
                                    grid_size=(14, 16),
                                    skip_SNR_above_hundred=True)



                # Grid of plots of x=SNR, y=knoise/kmax, across mach numbers
                print("----------------- FIGURE C -------------------------")
                plot_flows_grid(turb_data[res],
                                    res=res,
                                    flow_types=['s', 'c'],
                                    machs=mach_nums,
                                    grid_size=(9, 4),
                                    fit_params=fit_params)




                # Show the params of figure C varying.
                print("----------------- FIGURE D-------------------------")
                plot_fit_params_grid(turb_data[res],
                                    flow_types=['s', 'm', 'c'],
                                    machs=mach_nums,
                                    fit_params=fit_params)

                # figure C but over zetas
                print("----------------- FIGURE F -----------------------")
                plot_machs_grid(turb_data[res],
                                res=res,
                                machs=[0.1, 10],
                                grid_size=(9, 4),
                                fit_params=fit_params)

                # figure D but over zetas
                print("----------------- FIGURE G -------------------------")
                plot_fit_params_grid_zeta(turb_data[res],
                                            machs=mach_nums,
                                            fit_params=fit_params)



                print("----------------- FIGURE H -------------------------")
                plot_fig_H(turb_data[res],
                            machs=mach_nums,
                            fit_params=fit_params)


                print("----------------- FIGURE I -------------------------")
                plot_fig_I(turb_data[res],
                            machs=mach_nums,
                            fit_params=fit_params)


            if supersonic_mach_range:
                print(f"----------------- PLOTTING SUPERSONIC MACH RANGE FIGS ---------------------",
                        color='yellow')
                    # Grid of power spectra
                print("----------------- FIGURE B -------------------------")
                plot_spectra_grid(turb_data[res],
                                    res=res,
                                    flow_types=['s', 'm', 'c'],
                                    machs=[1, 2, 5, 10],
                                    grid_size=(14, 16),
                                    skip_SNR_above_hundred=True)



                # Grid of plots of x=SNR, y=knoise/kmax, across mach numbers
                print("----------------- FIGURE C -------------------------")
                plot_flows_grid(turb_data[res],
                                    res=res,
                                    flow_types=['s', 'c'],
                                    machs=supersonic_mach_nums,
                                    grid_size=(9, 4),
                                    fit_params=fit_params)




                # Show the params of figure C varying.
                print("----------------- FIGURE D-------------------------")
                plot_fit_params_grid(turb_data[res],
                                    flow_types=['s', 'm', 'c'],
                                    machs=supersonic_mach_nums,
                                    fit_params=fit_params)

                # figure C but over zetas
                print("----------------- FIGURE F -----------------------")
                plot_machs_grid(turb_data[res],
                                res=res,
                                machs=[1,10],
                                grid_size=(9, 4),
                                fit_params=fit_params)

                # figure D but over zetas
                print("----------------- FIGURE G -------------------------")
                plot_fit_params_grid_zeta(turb_data[res],
                                            machs=supersonic_mach_nums,
                                            fit_params=fit_params)



                print("----------------- FIGURE H -------------------------")
                plot_fig_H(turb_data[res],
                            machs=supersonic_mach_nums,
                            fit_params=fit_params)


                print("----------------- FIGURE I -------------------------")
                plot_fig_I(turb_data[res],
                            machs=supersonic_mach_nums,
                            fit_params=fit_params)


    print("End.")

