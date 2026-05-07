#!/usr/bin/env python

from globals_turb_data import *

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
def plot_noisy_map_examples(turb_data, void_size):
    # plot example maps of noisy data
    cfp.load_plot_style(fontsize=20)
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
        po.plt().text(0.5, 0.88, r"$\mathrm{SNR} = "+str(SNR)+r"$", ha='center', va='center',
                        transform=po.plt().gca().transAxes, bbox=dict(boxstyle="round", fc="white", alpha=0.9))
        cfp.show_or_save_plot(save=save)
    # plot colorbar
    cfp.plot_colorbar(vmin=vmin, vmax=vmax, cmap=cmap, label=r'Column density $(\mathrm{g\,cm^{-2}})$',
                        panels=3, orientation="horizontal", save=figpath+'noisy_data_map_cbar_'+void_size+'.pdf')
    load_default_plotstyle()

# ===================================================================
def plot_sigrho_Brunt_vs_SNR(turb_data, void_size):
    # figure 3 (Brunt sigma_rho as function of SNR without noise correction)
    xlim = [5e-3, 2e2]
    for isnr, SNR in enumerate(SNR_list):
        x = SNR
        y = turb_data['stats'][isnr].avg.Brunt_sigma_3D
        yerr = [y-turb_data['stats'][isnr].p16.Brunt_sigma_3D, turb_data['stats'][isnr].p84.Brunt_sigma_3D-y]
        cfp.plot(x=x, y=y, yerr=np.array([yerr]).T, color=void_colors[void_size],
                    marker=void_markers[void_size], alpha=0.9)
    # draw true sigma_rho
    y = turb_data['stats'][0].avg.sigma_3D
    yerr = [y-turb_data['stats'][0].p16.sigma_3D, turb_data['stats'][0].p84.sigma_3D-y]
    cfp.plot(x=xlim, y=[y,y], yerr=np.array([yerr,yerr]).T, shaded_err=True,
                linestyle=(0, (2, 2)), color=void_colors[void_size], label=r'True $\sigma_{\rho/\rho_0}$',
                legend_formatter=cfp.legend_formatter().add(pos=(0.63,0.89), length=1.35))
    # generate final plot and save
    cfp.plot(xlabel=r"$\mathrm{SNR}$", ylabel=r"$\sigma_{\rho/\rho_0,\mathrm{Brunt}}$",
                xlog=True, ylog=True, xlim=xlim, ylim=[1e-1,2e2],
                save=figpath+'Brunt_sigrho_vs_SNR_'+void_size+'.pdf')

# ===================================================================
def plot_spectra(turb_data, void_size, noise_exp=1.0, res=1024, res_resample=1.0):
    xlim = [1e-2, 2]
    ylim = [1e-14, 2e-2]
    kmax = turb_data['stats'][0].avg.k.max()
    x = turb_data['stats'][0].avg.k[1:]/kmax
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
        norm = turb_data['stats'][isnr].avg.PS_S[0]
        # draw the noise spectra
        y = turb_data['stats'][isnr].avg.PS_N[1:]
        stop()
        if SNR == 1: lf = cfp.legend_formatter().add(pos=(legend_xpos,legend_ypos-0.07), length=legend_length)
        else: lf = None
        cfp.plot(x=x, y=y, color='grey', linestyle=(0, (2, 2)), linewidth=1.0, label=labels[0], legend_formatter=lf)
        cfp.plot(x=1.1, y=10**(np.log10(y[-1])+pos_correct), text=str(SNR))
        # draw the noisy data spectra
        y = turb_data['stats'][isnr].avg.PS_SpN[1:]
        if SNR == 1: lf = cfp.legend_formatter().add(pos=(legend_xpos,legend_ypos), length=legend_length)
        else: lf = None
        po = cfp.plot(x=x, y=y, color=void_colors[void_size], label=labels[1], legend_formatter=lf)
        # draw k_noise
        k_noise = turb_data['stats'][isnr].avg.k_noise/kmax
        k_noise_min = turb_data['stats'][isnr].p16.k_noise/kmax
        k_noise_max = turb_data['stats'][isnr].p84.k_noise/kmax
        if SNR > 1e-2 and k_noise_min >= 2.5e-2 and k_noise_max <= 0.98: # exclude drawing very high SNR cases, as k_noise is effectively 1
            y_cross = turb_data['stats'][isnr].avg.PS_S_at_k_noise # intersection with signal spectrum
            po.plt().vlines(k_noise, ylim[0], y_cross, color=void_colors[void_size], linestyle=(0, (1, 1)))
            po.plt().fill_between([k_noise_min, k_noise_max], ylim[0], y_cross,
                                    color=void_colors[void_size], alpha=0.2, edgecolor='none')
            # fake line to get legend formatter to work
            if SNR == 1: lf = cfp.legend_formatter().add(pos=(legend_xpos+0.25,legend_ypos), length=legend_length)
            else: lf = None
            cfp.plot(x=0, y=1, color=void_colors[void_size], linestyle=(0, (1, 1)),
                        label=labels[2], legend_formatter=lf)
    # plot signal spectrum
    y = turb_data['stats'][0].avg.PS_S[1:]
    yerr = [y-turb_data['stats'][0].p16.PS_S[1:], turb_data['stats'][0].p84.PS_S[1:]-y]
    cfp.plot(x=x, y=y, yerr=yerr, shaded_err=True, color='red', label='S',
                legend_formatter=cfp.legend_formatter().add(pos=(legend_xpos,legend_ypos-0.14), length=legend_length))
    # generate final plot and save
    cfp.plot(x=0.886, y=0.92, text='SNR', normalised_coords=True)
    noise_str = ""
    if noise_exp != 1.0: noise_str = "_noise_exp_"+str(noise_exp)
    res_str = ""
    if res != 1024: res_str = "_res_"+str(res)
    res_resample_str = ""
    if res_resample != 1.0: res_resample_str = "_res_resample_"+str(res_resample)
    save = figpath+'PS_'+void_size+noise_str+res_str+res_resample_str+'.pdf'
    cfp.plot(xlabel=r"$k/k_\mathrm{max}$", ylabel=r"$P$", xlog=True, ylog=True, xlim=xlim, ylim=ylim, save=save)

# ===================================================================
def plot_k_noise_vs_SNR(turb_data):
    # figure 5 (k_noise vs SNR)
    for ivs, void_size in enumerate(void_sizes):
        kmax = turb_data[void_size]['stats'][0].avg.k.max()
        for isnr, SNR in enumerate(SNR_list):
            x = SNR
            y = turb_data[void_size]['stats'][isnr].avg.k_noise/kmax
            yerr = [y-turb_data[void_size]['stats'][isnr].p16.k_noise/kmax,
                    turb_data[void_size]['stats'][isnr].p84.k_noise/kmax-y]
            if SNR == 1:
                label = void_size_labels[void_size]
                lf = cfp.legend_formatter().add(pos=(0.02,0.89-ivs*0.09))
            else:
                label = None
                lf = None
            cfp.plot(x=x, y=y, yerr=np.array([yerr]).T, color=void_colors[void_size], linestyle='none',
                        marker=void_markers[void_size], alpha=0.8, label=label, legend_formatter=lf)
    # generate final plot and save
    cfp.plot(xlabel=r"$\mathrm{SNR}$", ylabel=r"$k_\mathrm{noise}/k_\mathrm{max}$",
                xlog=True, ylog=True, xlim=[5e-2,2e2], ylim=[2e-2,2e0], save=figpath+'k_noise_vs_SNR.pdf')

# ===================================================================
def plot_k_noise_vs_SNR_res(turb_data):
    # similar to figure 5 (k_noise vs SNR), but for fixed 50mu voids and resolutions 512, 1024, 2048
    resolutions = [512, 1024, 2048]
    res_colors = ['grey', 'blue', 'magenta']
    res_markers = ['v', 'o', '^']
    for ires, res in enumerate(resolutions):
        kmax = turb_data[str(res)]['stats'][0].avg.k.max()
        for isnr, SNR in enumerate(SNR_list):
            x = SNR
            y = turb_data[str(res)]['stats'][isnr].avg.k_noise/kmax
            yerr = [y-turb_data[str(res)]['stats'][isnr].p16.k_noise/kmax,
                    turb_data[str(res)]['stats'][isnr].p84.k_noise/kmax-y]
            if SNR == 1:
                label = r"$"+str(res)+r"$"
                lf = cfp.legend_formatter().add(pos=(0.02,0.89-ires*0.09))
            else:
                label = None
                lf = None
            cfp.plot(x=x, y=y, yerr=np.array([yerr]).T, color=res_colors[ires], linestyle='none',
                        marker=res_markers[ires], alpha=0.8, label=label, legend_formatter=lf)
    # generate final plot and save
    cfp.plot(xlabel=r"$\mathrm{SNR}$", ylabel=r"$k_\mathrm{noise}/k_\mathrm{max}$",
                xlog=True, ylog=True, xlim=[5e-2,2e2], ylim=[2e-2,2e0], save=figpath+'k_noise_vs_SNR_res.pdf')

# ===================================================================
def plot_k_noise_vs_SNR_res_resample(turb_data):
    # similar to figure 5 (k_noise vs SNR), but for fixed 50mu voids and resolution resampling factors 0.5, 1.0, 2.0
    res_resample = [0.5, 1.0, 2.0]
    res_colors = ['grey', 'blue', 'magenta']
    res_markers = ['v', 'o', '^']
    for ires, res in enumerate(res_resample):
        kmax = turb_data[str(res)]['stats'][0].avg.k.max()
        for isnr, SNR in enumerate(SNR_list):
            x = SNR
            y = turb_data[str(res)]['stats'][isnr].avg.k_noise/kmax
            yerr = [y-turb_data[str(res)]['stats'][isnr].p16.k_noise/kmax,
                    turb_data[str(res)]['stats'][isnr].p84.k_noise/kmax-y]
            if SNR == 1:
                label = r"Resampling factor $"+str(res)+r"$"
                lf = cfp.legend_formatter().add(pos=(0.40,0.45-ires*0.09))
            else:
                label = None
                lf = None
            cfp.plot(x=x, y=y, yerr=np.array([yerr]).T, color=res_colors[ires], linestyle='none',
                        marker=res_markers[ires], alpha=0.8, label=label, legend_formatter=lf)
    # generate final plot and save
    cfp.plot(xlabel=r"$\mathrm{SNR}$", ylabel=r"$k_\mathrm{noise}/k_\mathrm{max}$",
                xlog=True, ylog=True, xlim=[5e-2,2e2], ylim=[2e-2,2e0], save=figpath+'k_noise_vs_SNR_res_resample.pdf')

# ===================================================================
def plot_R_vs_SNR(turb_data):
    # figure 6 (R vs SNR, and R_denoised vs SNR)
    xlim = [5e-2, 2e2]
    plot_types = ["R_to_half", "R_to_half_denoised"]
    for ipt, pt in enumerate(plot_types):
        if pt == "R_to_half":
            ylim = [0.09, 0.45]
            xlabel = ""
            ylabel = r"$R^{1/2}$"
            axes_format=["", None]
        if pt == "R_to_half_denoised":
            ylim = [0.3, 0.66]
            xlabel = r"$\mathrm{SNR}$"
            ylabel = r"$R^{1/2}_\mathrm{denoised}$"
            axes_format=[None, None]
        for ivs, void_size in enumerate(void_sizes):
            for isnr, SNR in enumerate(SNR_list):
                x = SNR
                y = getattr(turb_data[void_size]['stats'][isnr].avg, pt)
                yerr = [y-getattr(turb_data[void_size]['stats'][isnr].p16, pt),
                        getattr(turb_data[void_size]['stats'][isnr].p84, pt)-y]
                # make legend
                if pt == "R_to_half_denoised" and SNR == 1:
                    label = void_size_labels[void_size]
                    lf = cfp.legend_formatter().add(pos=(0.49,0.90-ivs*0.09))
                else:
                    label = None
                    lf = None
                cfp.plot(x=x, y=y, yerr=np.array([yerr]).T, color=void_colors[void_size], linestyle='none',
                            marker=void_markers[void_size], alpha=0.8, label=label, legend_formatter=lf)
            # plot R^(1/2) at SNR->infinity
            y = getattr(turb_data[void_size]['stats'][0].avg, pt)
            yerr = [y-getattr(turb_data[void_size]['stats'][0].p16, pt),
                    getattr(turb_data[void_size]['stats'][0].p84, pt)-y]
            cfp.plot(x=xlim, y=[y,y], yerr=np.array([yerr,yerr]).T, shaded_err=True,
                        color=void_colors[void_size], linewidth=1.0, linestyle='dashed', alpha=0.8)
        # generate final plot and save
        if pt == "R_to_half_denoised": # fake line for R^(1/2) at SNR->infinity legend
            cfp.plot(x=[0.0,0.1], y=[0.0,0.0], color='grey', linewidth=1.0, linestyle='dashed', alpha=0.8,
                        label=r"$R^{1/2}(\mathrm{SNR}\to\infty)$",
                        legend_formatter=cfp.legend_formatter().add(pos=(0.49,0.90-0.29), length=1.4))
        cfp.plot(xlabel=xlabel, ylabel=ylabel, xlog=True, axes_format=axes_format,
                    xlim=xlim, ylim=ylim, save=figpath+pt+'_vs_SNR.pdf')

# ===================================================================
def plot_sigrho_Brunt_denoised_and_noise_corrections_vs_SNR(turb_data):
    # figure 7 (denoised Brunt sigma_rho and noise correctoions as function of SNR)
    xlim = [5e-2, 2e2]
    plot_types = ["Brunt_sigrho_denoised", "C1_C2", "C"]
    for ipt, pt in enumerate(plot_types):
        if pt == "Brunt_sigrho_denoised":
            ylim = [0.0, 1.0]
            ylog = False
            xlabel = ""
            ylabel = r"$\sigma_{\rho/\rho_0,\mathrm{Brunt,denoised}}$"
            axes_format = ["", None]
        if pt == "C1_C2":
            ylim = [5e-3, 5.0]
            ylog = True
            xlabel = ""
            ylabel = r"Noise correction parts $C_1$, $C_2$"
            axes_format = ["", None]
        if pt == "C":
            ylim = [1e-2, 2.0]
            ylog = True
            xlabel = r"$\mathrm{SNR}$"
            ylabel = r"Combined noise correction $C$"
            axes_format = [None, None]
        # --- gather data for each void size and SNR ---
        x = []
        sigrho_Brunt_denoised = []
        sigrho_Brunt_denoised_err = []
        sigrho_Brunt_denoised_std = []
        sigrho_Brunt = []
        sigrho_Brunt_err = []
        sigrho_Brunt_std = []
        sigrho_Brunt_denoised_at_SNR_inf = []
        sigrho_Brunt_denoised_at_SNR_inf_err = []
        sigrho_Brunt_denoised_at_SNR_inf_std = []
        C1 = []
        C1_std = []
        C2 = []
        C2_std = []
        C = []
        C_std = []
        for ivs, void_size in enumerate(void_sizes):
            # get sigma_rho Brunt denoised at SNR -> infinity
            dat = turb_data[void_size]['stats'][0].avg.Brunt_sigma_3D_denoised
            sigrho_Brunt_denoised_at_SNR_inf.append(dat)
            err = [dat-turb_data[void_size]['stats'][0].p16.Brunt_sigma_3D_denoised,
                    turb_data[void_size]['stats'][0].p84.Brunt_sigma_3D_denoised-dat]
            sigrho_Brunt_denoised_at_SNR_inf_err.append(err)
            ind_snr = 0
            for isnr, SNR in enumerate(SNR_list):
                if isnr < 2: continue # skip SNR -> infinity and SNR = 0.01
                x.append(SNR)
                # get sigma_rho Brunt denoised and its error
                dat = turb_data[void_size]['stats'][isnr].avg.Brunt_sigma_3D_denoised
                sigrho_Brunt_denoised.append(dat)
                err = [dat-turb_data[void_size]['stats'][isnr].p16.Brunt_sigma_3D_denoised,
                        turb_data[void_size]['stats'][isnr].p84.Brunt_sigma_3D_denoised-dat]
                sigrho_Brunt_denoised_err.append(err)
                # get sigma_rho Brunt and its error
                dat = turb_data[void_size]['stats'][isnr].avg.Brunt_sigma_3D
                sigrho_Brunt.append(dat)
                err = [dat-turb_data[void_size]['stats'][isnr].p16.Brunt_sigma_3D,
                        turb_data[void_size]['stats'][isnr].p84.Brunt_sigma_3D-dat]
                sigrho_Brunt_err.append(err)
                if pt == "C1_C2":
                    # C1 and error propagation
                    dat = sigrho_Brunt_denoised_at_SNR_inf[ivs] / sigrho_Brunt_denoised[ivs*6+ind_snr]
                    C1.append(dat)
                    std = np.mean(sigrho_Brunt_denoised_err[ivs*6+ind_snr])
                    sigrho_Brunt_denoised_std.append(std)
                    std_at_SNR_inf = np.mean(sigrho_Brunt_denoised_at_SNR_inf_err[ivs])
                    sigrho_Brunt_denoised_at_SNR_inf_std.append(std_at_SNR_inf)
                    C1std = dat * np.sqrt((std/sigrho_Brunt_denoised[ivs*6+ind_snr])**2 +
                                          (std_at_SNR_inf/sigrho_Brunt_denoised_at_SNR_inf[ivs])**2)
                    C1_std.append(C1std)
                    # C2 and error propagation
                    dat = sigrho_Brunt_denoised[ivs*6+ind_snr] / sigrho_Brunt[ivs*6+ind_snr]
                    C2.append(dat)
                    std = np.mean(sigrho_Brunt_denoised_err[ivs*6+ind_snr])
                    std_sigrho_Brunt = np.mean(sigrho_Brunt_err[ivs*6+ind_snr])
                    sigrho_Brunt_std.append(std_sigrho_Brunt)
                    C2std = dat * np.sqrt((std/sigrho_Brunt_denoised[ivs*6+ind_snr])**2 +
                                          (std_sigrho_Brunt/sigrho_Brunt[ivs*6+ind_snr])**2)
                    C2_std.append(C2std)
                if pt == "C":
                    dat = sigrho_Brunt_denoised_at_SNR_inf[ivs] / sigrho_Brunt[ivs*6+ind_snr]
                    C.append(dat)
                    std_at_SNR_inf = np.mean(sigrho_Brunt_denoised_at_SNR_inf_err[ivs])
                    std_sigrho_Brunt = np.mean(sigrho_Brunt_err[ivs*6+ind_snr])
                    Cstd = dat * np.sqrt((std_at_SNR_inf/sigrho_Brunt_denoised_at_SNR_inf[ivs])**2 +
                                         (std_sigrho_Brunt/sigrho_Brunt[ivs*6+ind_snr])**2)
                    C_std.append(Cstd)
                ind_snr += 1 # increase SNR index
        x = np.array(x)
        sigrho_Brunt_denoised = np.array(sigrho_Brunt_denoised)
        sigrho_Brunt_denoised_err = np.array(sigrho_Brunt_denoised_err)
        sigrho_Brunt_denoised_std = np.array(sigrho_Brunt_denoised_std)
        sigrho_Brunt = np.array(sigrho_Brunt)
        sigrho_Brunt_err = np.array(sigrho_Brunt_err)
        sigrho_Brunt_std = np.array(sigrho_Brunt_std)
        sigrho_Brunt_denoised_at_SNR_inf = np.array(sigrho_Brunt_denoised_at_SNR_inf)
        sigrho_Brunt_denoised_at_SNR_inf_err = np.array(sigrho_Brunt_denoised_at_SNR_inf_err)
        sigrho_Brunt_denoised_at_SNR_inf_std = np.array(sigrho_Brunt_denoised_at_SNR_inf_std)
        C1 = np.array(C1)
        C1_std = np.array(C1_std)
        C2 = np.array(C2)
        C2_std = np.array(C2_std)
        C = np.array(C)
        C_std = np.array(C_std)
        # --- end of getting data ---
        # plot data points
        for ivs, void_size in enumerate(void_sizes):
            if pt == "Brunt_sigrho_denoised":
                y = sigrho_Brunt_denoised[ivs*6:(ivs+1)*6]
                yerr = sigrho_Brunt_denoised_err[ivs*6:(ivs+1)*6].T
                legpos = [0.02, 0.90]
            if pt == "C1_C2":
                y = C1[ivs*6:(ivs+1)*6]
                yerr = np.array([C1_std[ivs*6:(ivs+1)*6], C1_std[ivs*6:(ivs+1)*6]])
                legpos = [0.58, 0.5]
            if pt == "C":
                y = C[ivs*6:(ivs+1)*6]
                yerr = np.array([C_std[ivs*6:(ivs+1)*6], C_std[ivs*6:(ivs+1)*6]])
                legpos = [0.58, 0.5]
            # plot data points
            cfp.plot(x=x[ivs*6:(ivs+1)*6], y=y, yerr=yerr, color=void_colors[void_size], linestyle='none',
                        marker=void_markers[void_size], alpha=0.8, label=void_size_labels[void_size],
                        legend_formatter=cfp.legend_formatter().add(pos=(legpos[0],legpos[1]-ivs*0.08)))
            # add C2 data points
            if pt == "C1_C2":
                y = C2[ivs*6:(ivs+1)*6]
                yerr = np.array([C2_std[ivs*6:(ivs+1)*6], C2_std[ivs*6:(ivs+1)*6]])
                cfp.plot(x=x[ivs*6:(ivs+1)*6], y=y, yerr=yerr, color=void_colors[void_size], linestyle='none',
                            marker=void_markers[void_size], alpha=0.8)
            # plot additional lines on Brunt_sigrho_denoised panel
            if pt == "Brunt_sigrho_denoised":
                # draw sigma_rho Brunt denoised at SNR -> infinity
                cfp.plot(x=xlim, y=[sigrho_Brunt_denoised_at_SNR_inf[ivs],sigrho_Brunt_denoised_at_SNR_inf[ivs]],
                            linewidth=1.0, linestyle=(0, (2, 2)), alpha=0.8, color=void_colors[void_size])
        # --- Fit functions for C1, C2, and C ---
        def C1_func(SNR, a, b):
            y = (SNR**a + b) / SNR**a
            return np.log(y)
        def C2_func(SNR, c, d):
            y = SNR**c / (SNR**c + d)
            return np.log(y)
        def C_func(SNR, a, b, c, d):
            y = C1_func(SNR,a,b) + C2_func(SNR,c,d)
            return y
        def get_fit_err_range(func, fitdat):
            from itertools import product
            xfit = cfp.logspace(xlim[0], xlim[1], 200)
            yfit = np.exp(func(xfit, *fitdat.popt)) # best-fit line
            # determine min/max fit range by looping over all parameters and error directions
            npar = len(fitdat.popt)
            yfit_tmp = []
            # loops over all combinations like (0,0), (0,1), (1,0), (1,1), ...
            for err_idx in product([0, 1], repeat=npar):
                pars = [fitdat.popt[i] + fitdat.perr[i][err_idx[i]] for i in range(npar)]
                yfit_tmp.append(func(xfit, *pars))
            yfit_tmp = np.exp(np.array(yfit_tmp))
            yfit_err = np.array([yfit-np.min(yfit_tmp,axis=0), np.max(yfit_tmp,axis=0)-yfit])
            return xfit, yfit, yfit_err
        # ---------------------------------------
        # fits to C1 and C2 on all voids combined
        if pt == "C1_C2":
            # C1
            y = np.log(C1)
            yerr = C1_std / C1
            params = {'a':[0.30,0.75,1.20], 'b':[0.01,0.10,0.30]}
            fitdat_C1 = cfp.fit(C1_func, x, y, yerr=yerr, params=params)
            xfit, yfit, yfit_err = get_fit_err_range(C1_func, fitdat_C1)
            cfp.plot(x=xfit, y=yfit, yerr=yfit_err, shaded_err=True, color='grey', linestyle='dotted', linewidth=1.0,
                        alpha=0.8, label=r"$C_1(\mathrm{SNR})$ fit",
                        legend_formatter=cfp.legend_formatter().add(pos=(legpos[0],legpos[1]-(ivs+1)*0.08-0.02), length=1.4))
            # C2
            y = np.log(C2)
            yerr = C2_std / C2
            params = {'c':[1.0,1.5,2.0], 'd':[2.0,3.0,4.0]}
            fitdat_C2 = cfp.fit(C2_func, x, y, yerr=yerr, params=params)
            xfit, yfit, yfit_err = get_fit_err_range(C2_func, fitdat_C2)
            cfp.plot(x=xfit, y=yfit, yerr=yfit_err, shaded_err=True, color='grey', linestyle='dashed', linewidth=1.0,
                        alpha=0.8, label=r"$C_2(\mathrm{SNR})$ fit",
                        legend_formatter=cfp.legend_formatter().add(pos=(legpos[0],legpos[1]-(ivs+2)*0.08-0.03), length=1.4))
        if pt == "C":
            # C
            fitdat_C = SimpleNamespace() # create fitdat_C obj and fill popt and perr for C
            fitdat_C.popt = [*fitdat_C1.popt, *fitdat_C2.popt]
            fitdat_C.perr = [*fitdat_C1.perr, *fitdat_C2.perr]
            xfit, yfit, yfit_err = get_fit_err_range(C_func, fitdat_C)
            cfp.plot(x=xfit, y=yfit, yerr=yfit_err, shaded_err=True, color='red', linestyle='solid', linewidth=1.0,
                        alpha=0.8, label=r"$C(\mathrm{SNR})$",
                        legend_formatter=cfp.legend_formatter().add(pos=(legpos[0],legpos[1]-(ivs+1)*0.08-0.02), length=1.4))
        # generate final plot and save
        if pt == "Brunt_sigrho_denoised": # fake line for SNR->infinity legend
            cfp.plot(x=[0.0,0.1], y=[0.0,0.0], color='grey', linewidth=1.0, linestyle=(0, (2, 2)), alpha=0.8,
                        label=r"$\sigma_{\rho/\rho_0,\mathrm{Brunt,denoised}}(\mathrm{SNR}\to\infty)$",
                        legend_formatter=cfp.legend_formatter().add(pos=(0.02,0.90-0.25), length=1.4))
        cfp.plot(xlabel=xlabel, ylabel=ylabel,
                    xlog=True, ylog=ylog, xlim=xlim, ylim=ylim, axes_format=["", None],
                    save=figpath+pt+'_vs_SNR.pdf')

# ===================================================================
def plot_noise_exponent_influence(turb_data_noise_exp):
    plot_types = ["k_noise", "Brunt_sigrho_denoised"]
    for pt in plot_types:
        if pt == "k_noise":
            xlabel = ""
            ylabel = r"$k_\mathrm{noise}/k_\mathrm{max}$"
            ylim = [-0.1, 1.15]
            axes_format = ["", None]
        if pt == "Brunt_sigrho_denoised":
            xlabel = r"Noise exponent $\alpha$"
            ylabel = r"$\sigma_{\rho/\rho_0,\mathrm{Brunt,denoised}}$"
            ylim = [0.08, 0.39]
            axes_format = [None, None]
        SNR_colors = ['grey', 'black', 'blue', 'green', 'orange', 'red']
        i_snr = 0
        for isnr, SNR in enumerate(SNR_list):
            if SNR < 0.1 or SNR > 100: continue
            noise_exponents = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
            y = []
            yerr = []
            for noise_exp in noise_exponents:
                if pt == "k_noise":
                    kmax = turb_data_noise_exp[str(noise_exp)]['stats'][isnr].avg.k.max()
                    dat = turb_data_noise_exp[str(noise_exp)]['stats'][isnr].avg.k_noise/kmax
                    err = [dat-turb_data_noise_exp[str(noise_exp)]['stats'][isnr].p16.k_noise/kmax,
                            turb_data_noise_exp[str(noise_exp)]['stats'][isnr].p84.k_noise/kmax-dat]
                if pt == "Brunt_sigrho_denoised":
                    dat = turb_data_noise_exp[str(noise_exp)]['stats'][isnr].avg.Brunt_sigma_3D_denoised
                    err = [dat-turb_data_noise_exp[str(noise_exp)]['stats'][isnr].p16.Brunt_sigma_3D_denoised,
                            turb_data_noise_exp[str(noise_exp)]['stats'][isnr].p84.Brunt_sigma_3D_denoised-dat]
                y.append(dat)
                yerr.append(err)
            cfp.plot(x=noise_exponents, y=y, yerr=np.array(yerr).T, color=SNR_colors[i_snr],
                        linestyle='solid', marker='o', alpha=0.8)
            if pt == "k_noise": cfp.plot(x=2.05, y=y[-1]-0.03, color=SNR_colors[i_snr], text=str(SNR))
            i_snr += 1
        # generate final plot and save
        if pt == "k_noise": cfp.plot(x=0.886, y=0.92, text='SNR', normalised_coords=True)
        cfp.plot(xlabel=xlabel, ylabel=ylabel, axes_format=axes_format,
                    xlim=[-0.2,2.3], ylim=ylim, save=figpath+pt+'_vs_noise_exp.pdf')


# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    load_default_plotstyle()

    # load default data
    turb_data = {}
    for void_size in void_sizes:
        dillfile = get_filename(1024, void_size, dill=True)
        with open(dillfile, "rb") as f:
            turb_data[void_size] = dill.load(f)
        print("Read '"+dillfile+"'.", color='green')

    # load noise exponent data
    turb_data_noise_exp = {}
    for noise_exp in [0.0, 0.5, 1.0, 1.5, 2.0]:
        dillfile = get_filename(1024, '50mu', noise_exp=noise_exp, dill=True)
        with open(dillfile, "rb") as f:
            turb_data_noise_exp[str(noise_exp)] = dill.load(f)
        print("Read '"+dillfile+"'.", color='green')

    # load simulation resolution data
    turb_data_res = {}
    for res in [512, 1024, 2048]:
        dillfile = get_filename(res, '50mu', dill=True)
        with open(dillfile, "rb") as f:
            turb_data_res[str(res)] = dill.load(f)
        print("Read '"+dillfile+"'.", color='green')

    # load resolution resampling data
    turb_data_res_resample = {}
    for res_resample in [1.0, 0.5, 2.0]:
        dillfile = get_filename(1024, '50mu', res_resample=res_resample, dill=True)
        with open(dillfile, "rb") as f:
            turb_data_res_resample[str(res_resample)] = dill.load(f)
        print("Read '"+dillfile+"'.", color='green')

    # k_noise vs SNR and spectra for different resolution resampling factors (0.5, 1.0, 2.0)
    plot_k_noise_vs_SNR_res_resample(turb_data_res_resample)
    plot_spectra(turb_data_res_resample['0.5'], '50mu', res_resample=0.5)
    plot_spectra(turb_data_res_resample['1.0'], '50mu', res_resample=1.0)
    plot_spectra(turb_data_res_resample['2.0'], '50mu', res_resample=2.0)

    stop()

    # k_noise vs SNR and spectra for different simulation resolutions (512, 1024, 2048)
    plot_k_noise_vs_SNR_res(turb_data_res)
    plot_spectra(turb_data_res['512'], '50mu', res=512)
    plot_spectra(turb_data_res['1024'], '50mu', res=1024)
    plot_spectra(turb_data_res['2048'], '50mu', res=2048)

    # figure 2 (noisy map examples)
    plot_noisy_map_examples(turb_data['50mu'], '50mu')

    # figure 3 (Brunt sigma_rho as function of SNR without noise correction)
    plot_sigrho_Brunt_vs_SNR(turb_data['50mu'], '50mu')

    # figure 4 (spectra)
    plot_spectra(turb_data['50mu'], '50mu')

    # figure 5 (k_noise vs SNR)
    plot_k_noise_vs_SNR(turb_data)

    # figure 6 (R vs SNR)
    plot_R_vs_SNR(turb_data)

    # figure 7 (sigma_rho Brunt denoised and noise correction factors vs SNR)
    plot_sigrho_Brunt_denoised_and_noise_corrections_vs_SNR(turb_data)

    # Appendix (figure 8)
    plot_spectra(turb_data['12p5mu'], '12p5mu')
    plot_spectra(turb_data['100mu'], '100mu')

    # Appendix (figure 9)
    plot_spectra(turb_data_noise_exp['0.0'], '50mu', noise_exp=0.0)
    plot_spectra(turb_data_noise_exp['2.0'], '50mu', noise_exp=2.0)

    # Appendix (figure 10)
    plot_noise_exponent_influence(turb_data_noise_exp)
