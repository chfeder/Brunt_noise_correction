from startlaser import *
#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Christoph Federrath and Isabella Gerrard, 2021-2022

import argparse
import shutil
import numpy as np
from astropy.io import fits
from scipy import fftpack
#from cfpack import hdfio, stop, print
#impor cfpack as cfp

#Gaussian kernel weights
def kernel_weights(grid, centre, FWHM):
    y = grid[0]
    x = grid[1]
    Rsq = (x-centre[0])**2 + (y-centre[1])**2
    sigma = FWHM/2.355
    kernel = np.exp(-Rsq/(2*sigma**2))
    return kernel

# compute weighted standard deviation
def my_std(data, weights=None, mask=None):
    # copy input data into work array
    d = np.copy(data)
    # deal with weights
    if weights is None:
        w = np.ones(data.shape)
    else:
        w = np.copy(weights)
    # apply mask
    if mask is not None:
        d = np.copy(d[mask])
        w = w[mask]
    # if there are any points in data after masking that still have NaN, return NaN
    if np.sum(d) == np.nan:
        return np.nan
    else:
        avg = np.average(d,    weights=w)
        msq = np.average(d**2, weights=w)
    return np.sqrt(msq-avg**2)

# compute Brunt et al. (2010) sqrt(R) factor, i.e., ratio of 2D-to-3D density dispersion
def get_Brunt_sqrt_R(power_spectrum_map, k, exclude_k0=True):
    # deal with wave number (k) input shape to get spacing of wave number map
    if k.shape[0] > 1 and k.shape[1] > 1:
        dkx = k[0,0]-k[0,1]
        dky = k[0,0]-k[1,0]
    else:
        dkx = 1
        dky = 1
    dk  = np.hypot(dkx,dky)
    # copy power spectrum
    P = np.copy(power_spectrum_map)
    # compute total of 2D and 3D power, excluding k=0 (excluding the mean)
    if exclude_k0: P[P.shape[0]//2, P.shape[1]//2] = 0.0
    P2D = 2 * np.pi * np.sum(P  ) * dk
    P3D = 4 * np.pi * np.sum(P*k) * dk
    Brunt_sqrt_R  = np.sqrt(P2D / P3D) # sqrt of ratio of 2D-to-3D power
    return Brunt_sqrt_R

# compute the 2D power spectrum of the moment-0 map
def get_power_spectrum(mom0, zero_pad=False, weights=None, mirror=False, normalise=False, show=False):
    if weights is None: # if weight is not provided, we make an array of ones; i.e., no change to the data
        weights = np.ones(mom0.shape)
    # construct map to do power spectrum analysis on
    scaled_map = np.copy(mom0)
    x_inds = np.argwhere(~np.isnan(scaled_map))[:,0]
    y_inds = np.argwhere(~np.isnan(scaled_map))[:,1]
    xmin = np.min(x_inds)
    xmax = np.max(x_inds)
    ymin = np.min(y_inds)
    ymax = np.max(y_inds)
    # select minimum bounding box that does not contain NaNs and is scaled by the SNR
    scaled_map = (scaled_map * weights)[xmin:xmax, ymin:ymax]
    indices_nan = np.where(np.isnan(scaled_map))
    scaled_map[indices_nan] = 0.0 # set NaNs to zero
    # apply mirroring to ensure periodicity of the moment-0 map before Fourier transform
    if mirror:
        mirrored_upperhalf = np.hstack((np.flipud(scaled_map),np.fliplr(np.flipud(scaled_map))))
        mirrored_lowerhalf = np.hstack((scaled_map,np.fliplr(scaled_map)))
        mirrored_map = np.vstack((mirrored_lowerhalf,mirrored_upperhalf))
        scaled_map = mirrored_map
    # zero padding (this is now default)
    if zero_pad:
        dims = scaled_map.shape
        scaled_map = np.pad(scaled_map, ((dims[0],dims[0]),(dims[1],dims[1])), 'constant')
    # make sure the data is float type
    data = scaled_map.astype(float)
    # pre-processing following Federrath et al. (2016); i.e., define fluctuations map
    if normalise: data = data/data.mean() - 1.0
    # 2D Fourier transformation
    ft = fftpack.fft2(data) / (data.shape[0]*data.shape[1])
    # center the transform so k = (0,0) is in the center
    ft_c = fftpack.fftshift(ft)
    # compute the power spectrum map
    power_spectrum_map = np.abs(ft_c*np.conjugate(ft_c))
    # undoing mirroring after Fourier transform
    if mirror:
        nx = power_spectrum_map.shape[0]
        ny = power_spectrum_map.shape[1]
        power_spectrum_map = power_spectrum_map[nx%2::2, ny%2::2]
    # define k (wave number) grid
    nx = power_spectrum_map.shape[0]
    ny = power_spectrum_map.shape[1]
    kx = np.linspace(-(nx//2), nx//2+(nx%2-1), nx)
    ky = np.linspace(-(ny//2), ny//2+(ny%2-1), ny)
    kx_grid, ky_grid = np.meshgrid(kx,ky,indexing="ij")
    k = np.hypot(kx_grid,ky_grid)
    # trim LHS (in case of even input data dimension; i.e., k=0 is in the center)
    power_spectrum_map = power_spectrum_map[(-(nx%2)+1):,(-(ny%2)+1):]
    k = k[(-(nx%2)+1):,(-(ny%2)+1):]
    if show:
        if show != "Fourier image": cfp.plot_map(data, show=True)
        else: cfp.plot_map(power_spectrum_map, log=True, cmap_label="Fourier power", show=True)
    return power_spectrum_map, k

# fits a set of 2D points, zin(x,y), to a plane (weights is the SNR)
def plane_gradient_fit(zin, weights=None):
    zin = np.array(zin)
    if weights is not None: weights = np.array(weights)
    # find the number of x values and the number of y values of zin
    xn = np.shape(zin)[1]
    yn = np.shape(zin)[0]
    # create x and y arrays corresponding to coordinates of z
    xin = np.tile(np.array(range(xn), dtype=float), (yn,1))
    yin = np.tile(np.array(range(yn), dtype=float).reshape(yn,1), (1, xn)).reshape(yn, xn)
    # replace indices for which zin is NaN with NaNs for x and y
    xin[np.isnan(zin)] = np.nan
    yin[np.isnan(zin)] = np.nan
    # get the effective size of z
    z_size = zin.size
    num_nans = (np.isnan(zin)).sum()
    n = z_size - num_nans
    # if there is not enough data, return the average of the given data
    if n < 3:
        print("too few points to fit plane; n = ", n, "; returning mean of data...")
        return np.nanmean(zin)
    # shift x, y, z to the center of gravity frame
    x0 = float(np.nansum(xin)) / n
    y0 = float(np.nansum(yin)) / n
    z0 = float(np.nansum(zin)) / n
    x = xin - x0
    y = yin - y0
    z = zin - z0
    # see if have the weights
    if weights is not None:
        s = np.nansum(weights)
        sx = np.nansum(weights * x)
        sy = np.nansum(weights * y)
        sz = np.nansum(weights * z)
        sxx = np.nansum(weights * x * x)
        syy = np.nansum(weights * y * y)
        sxy = np.nansum(weights * x * y)
        sxz = np.nansum(weights * x * z)
        syz = np.nansum(weights * y * z)
    else:
        s = float(n)
        sx = np.nansum(x)
        sy = np.nansum(y)
        sz = np.nansum(z)
        sxx = np.nansum(x * x)
        syy = np.nansum(y * y)
        sxy = np.nansum(x * y)
        sxz = np.nansum(x * z)
        syz = np.nansum(y * z)
    D = s * (sxx * syy - sxy * sxy) + sx * (sxy * sy - sx * syy) + sy * (sx * sxy - sxx * sy)
    r0 = sz * (sxx * syy - sxy * sxy) + sx * (sxy * syz - sxz * syy) + sy * (sxz * sxy - sxx * syz)
    r1 = s * (sxz * syy - sxy * syz) + sz * (sxy * sy - sx * syy)   + sy * (sx * syz - sxz * sy)
    r2 = s * (sxx * syz - sxz * sxy) + sx * (sxz * sy - sx * syz)   + sz * (sx * sxy - sxx * sy)
    r = np.array([r0, r1, r2]) / D
    # shift x, y, z back to original frame
    r[0] = r[0] + z0 - r[1] * x0 - r[2] * y0
    plane_fit = r[0] + r[1] * xin + r[2] * yin
    return plane_fit

# fit gradient (plane) to log10(moment-0) map (using SNR as weights)
# return fitted plane and gradient-corrected moment-0
def get_corrected_moment0(mom0, weights=None):
    mom0_norm = mom0/np.nanmean(mom0)
    mom0_fit = 10**(plane_gradient_fit(np.log10(mom0_norm), weights=weights))
    mom0_corrected = mom0_norm / mom0_fit
    return mom0_norm, mom0_fit, mom0_corrected

# fit gradient (plane) to moment-1 map (using SNR as weights)
# return fitted plane and gradient-corrected moment-1
def get_corrected_moment1(mom1, weights=None):
    mom1_fit = plane_gradient_fit(mom1, weights=weights)
    mom1_corrected = mom1 - mom1_fit
    return mom1_fit, mom1_corrected

# calculate the turbulence driving parameter b = sigam_rho / Mach,
# using input moment maps (mom-0 is normalised internally, but mom-1 must be in km/s)
# and sound speed in km/s; SNR is the signal-to-noise-ratio map (same size as moment maps)
def turbulence_driving_analysis(mom0, mom1, snr_map=None, kernel=None, sound_speed=8.3, plot=False):
    # fit moment-0 with plane gradient
    mom0_norm, mom0_fit, mom0_corrected = get_corrected_moment0(mom0, weights=snr_map)
    # fit moment-1 with plane gradient
    mom1_fit, mom1_corrected = get_corrected_moment1(mom1, weights=snr_map)
    # get 2D power spectrum of moment-0
    power_spectrum_map, k = get_power_spectrum(mom0_corrected, weights=snr_map, mirror=False)
    # plot
    if plot:
        cfp.plot_map(mom0, save='mom0.pdf', log=True)
        cfp.plot_map(mom0_norm, save='mom0_norm.pdf', log=True)
        cfp.plot_map(mom0_fit, save='mom0_fit.pdf', log=True)
        cfp.plot_map(mom0_corrected, save='mom0_corrected.pdf', log=True)
        cfp.plot_map(mom1, save='mom1.pdf')
        cfp.plot_map(mom1_fit, save='mom1_fit.pdf')
        cfp.plot_map(mom1_corrected, save='mom1_corrected.pdf')
        cfp.plot_map(power_spectrum_map, save='power_spectrum_map.pdf', log=True)
    # apply Brunt et al. 2010 method to convert 2D (column) density dispersion to 3D (volume) density dispersion
    Brunt_sqrt_R = get_Brunt_sqrt_R(power_spectrum_map, k)
    # get 2D density dispersion
    sigma_rho_2D = my_std(mom0_corrected, weights=snr_map, mask=kernel)
    # get 3D density dispersion
    sigma_rho_3D = sigma_rho_2D / Brunt_sqrt_R
    # get 1D velocity dispersion
    sigma_vel_1D = my_std(mom1_corrected, weights=snr_map, mask=kernel)
    # get 3D velocity dispersion, using 1D-to-3D velocity dispersion correction factor of 3.3 for mom-1 from Stewart & Federrath (2022)
    sigma_vel_3D = 3.3 * sigma_vel_1D
    # get the turbulence driving parameter b = sigma_rho / Mach (all 3D)
    turbulence_driving_parameter_b = sigma_rho_3D / (sigma_vel_3D / sound_speed)
    return turbulence_driving_parameter_b

# generate turbulent field
def generate(args):
    if args.turb_file is None:
        if shutil.which("TurbGen") is None: print("Need TurbGen executable in path; go to 'https://www.mso.anu.edu.au/~chfeder/codes.html' for details.", error=True)
        args.turb_file = "turbulence_driving_analysis.h5" # default filename for this script
        print("Calling TurbGen to generate turbulent field in file '"+args.turb_file+"'...", highlight=1)
        arg_N = ""
        if args.ncells: arg_N = " -N "+str(args.ncells)+" "+str(args.ncells)+" "+str(args.ncells)
        arg_kmin = " -kmin "+str(args.kmin)
        arg_kmax = " -kmax "+str(args.kmax)
        args_spect_form = " -spect_form "+str(args.spect_form)
        args_power_law_exp = ""
        if args.spect_form == 2: args_power_law_exp = " -power_law_exp "+str(args.power_law_exp)
        args_angles_exp = " -angles_exp "+str(args.angles_exp)
        args_random_seed = " -random_seed "+str(args.random_seed)
        args_outputfile = " -o "+args.turb_file
        args_verbose = " -verbose "+str(args.verbose)
        cmd = "TurbGen -ndim 3"+arg_N+arg_kmin+arg_kmax+args_spect_form+args_power_law_exp+args_angles_exp
        cmd += args_random_seed+args_verbose+args_outputfile
        cfp.run_shell_command(cmd) # run TurbGen
    # read back
    dat = hdfio.read(args.turb_file, "turb_field_x") # note that mean(dat) = 0 and std(dat) = 1.0
    # now transform Gaussian x component to become a log-normal density field
    s = dat * args.sigma_s # scale by target sigma_s
    s = s - 0.5*s.std()**2 # first impose relation <s> = - sigma_s^2 / 2
    print("mean and std of s=ln(rho/<rho>) = ", s.mean(), s.std())
    pdf, s_bins = cfp.get_pdf(s) # PDF of s
    cfp.plot(x=s_bins, y=pdf, label='before correction')
    rho = np.exp(s) # define density field
    print("mean and std of rho/<rho> (before correction) = ", rho.mean(), rho.std())
    # apply small correction to make the mean exactly 0
    rho /= rho.mean()
    print("mean and std of rho/<rho> (after correction) = ", rho.mean(), rho.std())
    pdf, s_bins = cfp.get_pdf(np.log(rho)) # recompute PDF of s
    cfp.plot(x=s_bins, y=pdf, label='after correction', xlabel=r"$s$", ylabel=r"$\mathrm{PDF}(s)$", ylog=False, show=True)
    return rho

# compute column density along axis (dlos is the cell size)
def get_column_density_map(rho, dlos, axis=2):
    coldens = np.sum(rho, axis=axis) * dlos
    return coldens

# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Turbulence analysis.')

    parser.add_argument("-turb_file",
        help="File containing turbulent density field (generated by TurbGen; if not provided, generate via following options).")
    parser.add_argument("-ncells", type=int, default=128,
        help="Number of cells (default: %(default)s).")
    parser.add_argument('-kmin', type=float, default=2.0,
        help='Minimum wavenumber of generated field in units of 2pi/L[X]; (default: %(default)s).')
    parser.add_argument('-kmax', type=float, default=20.0,
        help='Maximum wavenumber of generated field in units of 2pi/L[X]; (default: %(default)s).')
    parser.add_argument('-spect_form', type=int, choices=[0, 1, 2], default=2,
        help='Spectral form: 0 (band/rectangle/constant), 1 (paraboloid), 2 (power law); (default: %(default)s).')
    parser.add_argument('-power_law_exp', type=float, default=-2.0,
        help="If spect_form 2: power-law exponent of the s=ln(rho/<rho>) spectrum " + \
                "(see e.g., Federrath et al. (2010), Fig. 15, for examples; (default: %(default)s).")
    parser.add_argument('-angles_exp', type=float, default=0.5,
        help='If spect_form 2: angles exponent for sparse sampling (e.g., 2.0: full sampling, 0.0: healpix-like sampling); (default: %(default)s).')
    parser.add_argument('-random_seed', type=int, default=140281,
        help='Random seed for turbulent field; (default: %(default)s).')

    parser.add_argument("-sigma_s", type=float, default=1.0,
        help="Target sigma_s (standard deviation of s=ln(rho/<rho>)) (default: %(default)s).")

    parser.add_argument("-ncubes", type=int, default=1,
        help="Number of replications of the turbulent cube along the LOS, to simulate an extended LOS depth of the medium (default: %(default)s).")

    parser.add_argument("-axis", type=int, default=2,
        help="Axis along which to integrate the volume density to obtain the column density (default: %(default)s).")

    parser.add_argument("-verbose", "--verbose", type=int, choices=[0, 1, 2], default=1,
        help="0 (no shell output), 1 (standard shell output), 2 (more shell output); (default: %(default)s)")

    args = parser.parse_args()

    # generate 3D density field (a cube)
    dens_cube = generate(args)

    # produce replicas of the parent cube along the LOS
    dens = np.concatenate([dens_cube]*args.ncubes, axis=args.axis)

    #dens = np.concatenate([dens_cube, np.flip(dens_cube, axis=1)], axis=args.axis)

    # show a slice (in the middle of the cube); prependicular to the LOS
    mid_index = np.s_[dens_cube.shape[0]//2]
    ind = [np.s_[:]]*3
    ind[args.axis-1] = mid_index # slice through middle of axis-1, so we see the replicated axis
    dens_slice = dens[tuple(ind)]
    cfp.plot_map(dens_slice, log=True, cmap_label="density slice", show=True)

    # volume density dispersion
    sigma_3D = np.std(dens)
    print("sigma_3D = ", sigma_3D)

    # compute cell size
    dlos = 1.0 / dens_cube.shape[0] # size of cube 1 x 1 x 1

    # compute column density
    coldens = get_column_density_map(dens, dlos, axis=args.axis)
    cfp.plot_map(coldens, log=True, cmap_label="column density", show=True)

    # normalise by the mean column density
    coldens /= coldens.mean()

    # column density dispersion
    sigma_2D = np.std(coldens)
    print("sigma_2D = ", sigma_2D)

    # get 2D power spectrum of column density map
    power_spectrum_map, k = get_power_spectrum(coldens, zero_pad=False, weights=None, mirror=False, normalise=False, show='Fourier image')

    # apply Brunt et al. (2010) method to convert 2D (column) density dispersion to 3D (volume) density dispersion
    Brunt_sqrt_R = get_Brunt_sqrt_R(power_spectrum_map, k)
    print("Brunt R^(1/2) = ", Brunt_sqrt_R)
    sigma_3D_from_Brunt = sigma_2D / Brunt_sqrt_R
    print("sigma_3D_from_Brunt = ", sigma_3D_from_Brunt)
