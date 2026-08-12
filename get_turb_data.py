#!/usr/bin/env python

from globals_turb_data import *

# ===================================================================
def get_turb_box(filename, box, res_resample=1.0):
    import yt
    yt.funcs.mylog.setLevel("ERROR")
    # box dimensions
    delta_a = 0.02
    edge_mtop = 0.35+delta_a
    edge_mbot = 0.35-delta_a
    # read FLASH file and extract box region via yt
    ds = yt.load(filename)
    if box == 'c_box':
       ledge = ds.arr([-delta_a, edge_mbot, -delta_a],'code_length')
       redge = ds.arr([+delta_a, edge_mtop, +delta_a],'code_length')
    if box == 'r_box':
        ledge = ds.arr([delta_a, edge_mbot, -delta_a],'code_length')
        redge = ds.arr([delta_a+2*delta_a, edge_mtop, delta_a],'code_length')
    if box == 'f_box':
        ledge = ds.arr([-delta_a, edge_mbot, delta_a],'code_length')
        redge = ds.arr([+delta_a, edge_mtop, delta_a+2*delta_a],'code_length')
    if box == 'l_box':
        ledge = ds.arr([-delta_a-2*delta_a, edge_mbot, -delta_a],'code_length')
        redge = ds.arr([-delta_a, edge_mtop, delta_a],'code_length')
    if box == 'b_box':
        ledge = ds.arr([-delta_a, edge_mbot, -delta_a-2*delta_a],'code_length')
        redge = ds.arr([+delta_a, edge_mtop, -delta_a],'code_length')
    dx = ds.index.get_smallest_dx()
    size = ((redge-ledge)/dx).astype('int')
    level = 0
    section = ds.covering_grid(level, ledge, size)
    dens = section['density'].value
    # requested resolution resampling by res_resample factor
    if res_resample != 1.0:
        new_shape = tuple(np.round(np.array(dens.shape)*res_resample).astype(int))
        print('Requested -- reshaping ', dens.shape, ' -> ', new_shape, warn=True)
        mean, std = dens.mean(), dens.std()
        dens = cfp.congrid(dens, new_shape)
        dens = (dens - dens.mean())/dens.std()*std + mean # correct to get exact previous mean and std
    # treat odd shape by interpolation
    if (np.array(dens.shape)%2)[0] != 0:
        new_shape = tuple(np.array(dens.shape)-1)
        print('Odd shape -- reshaping ', dens.shape, ' -> ', new_shape, warn=True)
        mean, std = dens.mean(), dens.std()
        dens = cfp.congrid(dens, new_shape)
        dens = (dens - dens.mean())/dens.std()*std + mean # correct to get exact previous mean and std
    return dens # density on uniform grid


# ===== the following applies in case we are running this in script mode =====
if __name__ == "__main__":

    # user changes here
    resolutions = [1024, 512, 2048] # simulation resolution (default first)
    void_sizes = ['50mu', '12p5mu', '100mu'] # foam void sizes (default first)
    noise_samples = 100 # convergence was checked also with 1000 noise samples, which yields practically the same
    print('Number of noise samples: '+str(noise_samples), color='lightred')

    timer = cfp.timer("MAIN")

    # loop over simulation resolutions
    for res in resolutions:

        # loop over foam void sizes
        for void_size in void_sizes:
            if res != 1024 and void_size != '50mu': continue # skip 512 and 2048 resolutions for 12.5 and 100 mu
            print("===== SIMULATION RESOLUTION AND VOID SIZE = "+str(res)+", "+void_size+" =====", color='lightred')

            # loop over noise exponents
            for noise_exp in [1.0, 0.0, 2.0, 0.5, 1.5]:

                # skip noise exponents other than 1 (white noise) for foam sizes 12.5mu and 100mu
                if noise_exp != 1.0 and (void_size != '50mu' or res != 1024): continue

                noise_exp_info = str(noise_exp)
                if noise_exp == 1.0: noise_exp_info += ' (white noise)'
                print("===== NOISE EXPONENT = "+noise_exp_info+" =====", color='cyan')

                # loop over resolution resampling cases
                for res_resample in res_resample_list:
                    if res_resample != 1.0 and res != 1024: continue

                    if res_resample != 1.0:
                        print("===== RESOLUTION RESAMPLING FACTOR = "+str(res_resample)+" =====", color='yellow')

                    # data collection per box and/or LOS
                    sigma_3D = []
                    data = [[] for _ in range(len(SNR_list))]

                    filename = get_filename(res, void_size)

                    # loop over the 5 boxes
                    for box in BOX_list:

                        dens = get_turb_box(filename, box, res_resample=res_resample)
                        print("=== Read "+box+" with dims "+str(dens.shape)+" from file "+filename, color='lightblue')

                        sigma_3D = (dens/dens.mean()).std()

                        # loop over the 3 projection axes (x,y,z)
                        for LOS in LOS_list:
                            print("computing projection along axis "+str(LOS), color='lightblue')
                            # five box locations
                            delta_a = 0.02
                            proj = dens.sum(axis=LOS) / dens.shape[LOS] * 2*delta_a

                            # loop over the 7 SNRs
                            for iSNR, SNR in enumerate(SNR_list):

                                ret = noise_iterations(noise_samples, proj, SNR, noise_exp=noise_exp)
                                ret.sigma_3D = np.array([sigma_3D]) # attach sigma_3D
                                data[iSNR].append(ret)

                    # compute statistics over noise samples, boxes and LOS'
                    stats = []
                    for isnr in range(len(SNR_list)):
                        avg, std, p16, p84 = get_stats_object(data[isnr])
                        stats.append(SimpleNamespace(avg=avg, std=std, p16=p16, p84=p84))

                    # package and save (dill) all data
                    turb_data = {"data": data, "stats": stats}
                    dillfile = get_filename(res, void_size, noise_exp=noise_exp,
                                            res_resample=res_resample, dill=True)
                    with open(dillfile, "wb") as f:
                        dill.dump(turb_data, f)
                    print("Wrote '"+dillfile+"'.", color='magenta')

    timer.report()