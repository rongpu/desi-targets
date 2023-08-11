from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp

sys.path.append(os.path.expanduser('~/git/desi-examples/imaging_systematics'))
from plot_healpix_map import plot_map


params = {'legend.fontsize': 'x-large',
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)

plt.rcParams['image.cmap'] = 'jet'

target_class = 'LRG'
# target_class = 'ELG_LOP'

if target_class=='BGS_BRIGHT':
    target_class = 'BGS_ANY'
    sub_class = 'BGS_BRIGHT'
elif target_class=='ELG_LOP':
    target_class = 'ELG'
    sub_class = 'ELG_LOP'
else:
    sub_class = target_class

min_nobs = 1
# maskbits = sorted([1, 13])
# maskbits = sorted([1, 12, 13])
# maskbits = sorted([1, 11, 12, 13])
# maskbits = sorted([1, 8, 9, 11, 12, 13])
# custom_mask_name = ''

maskbits = []
custom_mask_name = 'lrgmask_v1.1'
# custom_mask_name = 'elgmask_v1'

mask_str = ''.join([str(tmp) for tmp in maskbits])
if custom_mask_name!='':
    mask_str += '_' + custom_mask_name

randoms_counts_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/resolve/counts'
systematics_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/density_maps/1.1.1/resolve'
# target_densities_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/unofficial/density_maps'

top_plot_dir = '/global/cfs/cdirs/desicollab/users/rongpu/imaging_sys/systematics_maps/targets'

min_pix_frac = 0.2  # minimum fraction of pixel area to be used

xnames = ['FRACAREA', 'EBV', 'galdepth_gmag', 'galdepth_rmag', 'galdepth_zmag', 'psfdepth_gmag', 'psfdepth_rmag', 'psfdepth_zmag', 'psfdepth_w1mag', 'psfdepth_w2mag', 'galdepth_gmag_ebv', 'galdepth_rmag_ebv', 'galdepth_zmag_ebv', 'psfdepth_gmag_ebv', 'psfdepth_rmag_ebv', 'psfdepth_zmag_ebv', 'psfdepth_w1mag_ebv', 'psfdepth_w2mag_ebv', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 'NOBS_G', 'NOBS_R', 'NOBS_Z']
# xlabels = ['Heapix occupation fraction', 'GAIA stellar density [deg$^{-2}$]', 'log10(GAIA stellar density [deg$^{-2}$])', 'E(B-V)', 'g-band galaxy depth [mag]', 'r-band galaxy depth [mag]', 'z-band galaxy depth [mag]', 'g-band PSF depth [mag]', 'r-band PSF depth [mag]', 'z-band PSF depth [mag]', 'W1-band PSF depth [mag]', 'W2-band PSF depth [mag]', 'g-band galaxy depth - 3.214*E(B-V) [mag]', 'r-band galaxy depth - 2.165*E(B-V) [mag]', 'z-band galaxy depth - 1.211*E(B-V) [mag]', 'g-band PSF depth - 3.214*E(B-V) [mag]', 'r-band PSF depth - 2.165*E(B-V) [mag]', 'z-band PSF depth - 1.211*E(B-V) [mag]', 'W1-band PSF depth - 0.184*E(B-V) [mag]', 'W2-band PSF depth - 0.113*E(B-V) [mag]', 'g-band PSF size [arcsec]', 'r-band PSF size [arcsec]', 'z-band PSF size [arcsec]', 'NOBS_G', 'NOBS_R', 'NOBS_Z']
xlabels = ['Heapix occupation fraction', 'E(B-V)', 'g-band galaxy depth [mag]', 'r-band galaxy depth [mag]', 'z-band galaxy depth [mag]', 'g-band PSF depth [mag]', 'r-band PSF depth [mag]', 'z-band PSF depth [mag]', 'W1-band PSF depth [mag]', 'W2-band PSF depth [mag]', 'g-band galaxy depth - 3.214*E(B-V) [mag]', 'r-band galaxy depth - 2.165*E(B-V) [mag]', 'z-band galaxy depth - 1.211*E(B-V) [mag]', 'g-band PSF depth - 3.214*E(B-V) [mag]', 'r-band PSF depth - 2.165*E(B-V) [mag]', 'z-band PSF depth - 1.211*E(B-V) [mag]', 'W1-band PSF depth - 0.184*E(B-V) [mag]', 'W2-band PSF depth - 0.113*E(B-V) [mag]', 'g-band PSF size [arcsec]', 'r-band PSF size [arcsec]', 'z-band PSF size [arcsec]', 'NOBS_G', 'NOBS_R', 'NOBS_Z']

bin_params = {}
bin_params['FRACAREA'] = np.array([0., 1.])
bin_params['EBV'], bin_params['EBV_nbins'] = np.array([0., 0.15]), 30
bin_params['galdepth_gmag'], bin_params['galdepth_gmag_nbins'] = np.array([23.75, 25.0]), 50
bin_params['galdepth_rmag'], bin_params['galdepth_rmag_nbins'] = np.array([23.1, 24.8]), 50
bin_params['galdepth_zmag'], bin_params['galdepth_zmag_nbins'] = np.array([22.5, 23.65]), 50
bin_params['psfdepth_gmag'], bin_params['psfdepth_gmag_nbins'] = np.array([23.75+0.15, 25.0+0.35]), 50
bin_params['psfdepth_rmag'], bin_params['psfdepth_rmag_nbins'] = np.array([23.1+0.2, 24.8+0.4]), 50
bin_params['psfdepth_zmag'], bin_params['psfdepth_zmag_nbins'] = np.array([22.5+0.2, 23.65+0.4]), 50
bin_params['psfdepth_w1mag'], bin_params['psfdepth_w1mag_nbins'] = np.array([21.2, 22.2]), 50
bin_params['psfdepth_w2mag'], bin_params['psfdepth_w2mag_nbins'] = np.array([20.5, 21.8]), 50
bin_params['PSFSIZE_G'], bin_params['PSFSIZE_G_nbins'] = np.array([1.1, 2.5]), 30
bin_params['PSFSIZE_R'], bin_params['PSFSIZE_R_nbins'] = np.array([1.0, 2.3]), 30
bin_params['PSFSIZE_Z'], bin_params['PSFSIZE_Z_nbins'] = np.array([0.9, 1.9]), 30
bin_params['NOBS_G'], bin_params['NOBS_G_nbins'] = np.array([0, 8]), 8
bin_params['NOBS_R'], bin_params['NOBS_R_nbins'] = np.array([0, 8]), 8
bin_params['NOBS_Z'], bin_params['NOBS_Z_nbins'] = np.array([0, 8]), 8
bin_params['galdepth_gmag_ebv'], bin_params['galdepth_gmag_ebv_nbins'] = bin_params['galdepth_gmag'] - 0.1, 50
bin_params['galdepth_rmag_ebv'], bin_params['galdepth_rmag_ebv_nbins'] = bin_params['galdepth_rmag'] - 0.1, 50
bin_params['galdepth_zmag_ebv'], bin_params['galdepth_zmag_ebv_nbins'] = bin_params['galdepth_zmag'] - 0.05, 50
bin_params['psfdepth_w1mag_ebv'], bin_params['psfdepth_w1mag_ebv_nbins'] = bin_params['psfdepth_w1mag'], 50
bin_params['psfdepth_w2mag_ebv'], bin_params['psfdepth_w2mag_ebv_nbins'] = bin_params['psfdepth_w2mag'], 50
bin_params['psfdepth_gmag_ebv'], bin_params['psfdepth_gmag_ebv_nbins'] = bin_params['psfdepth_gmag'] - 0.1, 50
bin_params['psfdepth_rmag_ebv'], bin_params['psfdepth_rmag_ebv_nbins'] = bin_params['psfdepth_rmag'] - 0.1, 50
bin_params['psfdepth_zmag_ebv'], bin_params['psfdepth_zmag_ebv_nbins'] = bin_params['psfdepth_zmag'] - 0.05, 50
bin_params['psfdepth_w1mag_ebv'], bin_params['psfdepth_w1mag_ebv_nbins'] = bin_params['psfdepth_w1mag'], 50
bin_params['psfdepth_w2mag_ebv'], bin_params['psfdepth_w2mag_ebv_nbins'] = bin_params['psfdepth_w2mag'], 50

# for index in range(len(xnames)):
#     print(xnames[index], xlabels[index])

# nside = 64
for nside in [64, 128, 256]:

    npix = hp.nside2npix(nside)
    pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
    print(nside, 'Healpix size = {:.5f} sq deg'.format(pix_area))

    field = 'combined'

    plot_dir = os.path.join(top_plot_dir, '{}_minobs_{}_maskbits_{}'.format(field, min_nobs, mask_str))
    if not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)

    maps = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('north', nside, min_nobs, mask_str)))
    maps = maps[maps['n_randoms']>0]
    maps1 = Table.read(os.path.join(systematics_dir,'systematics_map_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(sub_class.lower(), 'north', nside, min_nobs, mask_str)))
    maps1.remove_columns(['RA', 'DEC'])
    maps = join(maps, maps1, join_type='inner', keys='HPXPIXEL')
    maps_north = maps.copy()

    maps = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('south', nside, min_nobs, mask_str)))
    maps = maps[maps['n_randoms']>0]
    maps1 = Table.read(os.path.join(systematics_dir, 'systematics_map_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(sub_class.lower(), 'south', nside, min_nobs, mask_str)))
    maps1.remove_columns(['RA', 'DEC'])
    maps = join(maps, maps1, join_type='inner', keys='HPXPIXEL')
    maps_south = maps.copy()

    mask = (maps_north['DEC']>32.375)
    maps_north = maps_north[mask]

    mask = ~np.in1d(maps_south['HPXPIXEL'], maps_north['HPXPIXEL'])
    maps = vstack([maps_north, maps_south[mask]])

    print(len(maps))

    area = np.sum(maps['FRACAREA'])*pix_area
    print('Area = {:.1f} sq deg'.format(area))

    mask = maps['FRACAREA']>min_pix_frac
    maps = maps[mask]

    # # Load stellar density map
    # stardens = np.load('/global/cfs/cdirs/desi/users/rongpu/useful/healpix_maps/pixweight-dr7.1-0.22.0_stardens_{}_ring.npy'.format(nside))
    # maps['stardens'] = stardens[maps['HPXPIXEL']]
    # maps['stardens_log'] = np.log10(maps['stardens'])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        maps['galdepth_gmag'] = -2.5*(np.log10((5/np.sqrt(maps['GALDEPTH_G'])))-9)
        maps['galdepth_rmag'] = -2.5*(np.log10((5/np.sqrt(maps['GALDEPTH_R'])))-9)
        maps['galdepth_zmag'] = -2.5*(np.log10((5/np.sqrt(maps['GALDEPTH_Z'])))-9)
        maps['psfdepth_gmag'] = -2.5*(np.log10((5/np.sqrt(maps['PSFDEPTH_G'])))-9)
        maps['psfdepth_rmag'] = -2.5*(np.log10((5/np.sqrt(maps['PSFDEPTH_R'])))-9)
        maps['psfdepth_zmag'] = -2.5*(np.log10((5/np.sqrt(maps['PSFDEPTH_Z'])))-9)
        maps['psfdepth_w1mag'] = -2.5*(np.log10((5/np.sqrt(maps['PSFDEPTH_W1'])))-9)
        maps['psfdepth_w2mag'] = -2.5*(np.log10((5/np.sqrt(maps['PSFDEPTH_W2'])))-9)
        maps['galdepth_gmag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps['GALDEPTH_G'])))-9) - 3.214*maps['EBV']
        maps['galdepth_rmag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps['GALDEPTH_R'])))-9) - 2.165*maps['EBV']
        maps['galdepth_zmag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps['GALDEPTH_Z'])))-9) - 1.211*maps['EBV']
        maps['psfdepth_gmag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps['PSFDEPTH_G'])))-9) - 3.214*maps['EBV']
        maps['psfdepth_rmag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps['PSFDEPTH_R'])))-9) - 2.165*maps['EBV']
        maps['psfdepth_zmag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps['PSFDEPTH_Z'])))-9) - 1.211*maps['EBV']
        maps['psfdepth_w1mag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps['PSFDEPTH_W1'])))-9) - 0.184*maps['EBV']
        maps['psfdepth_w2mag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps['PSFDEPTH_W2'])))-9) - 0.113*maps['EBV']

    for index, xname in enumerate(xnames):

        plot_path = os.path.join(plot_dir, 'map_{}_{}_{}.png'.format(sub_class.lower(), xname, nside))

        if os.path.isfile(plot_path):
            continue

        plot_map(nside, maps[xname], pix=maps['HPXPIXEL'],
                 vmin=bin_params[xname][0], vmax=bin_params[xname][1],
                 title='{} NSIDE={}'.format(xlabels[index], nside), save_path=plot_path, show=False)
