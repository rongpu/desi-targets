from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp

params = {'legend.fontsize': 'x-large',
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)

plt.rcParams['image.cmap'] = 'jet'

target_ver_str = '1.0.0'

min_nobs = 1
# maskbits_dict = {'LRG': [1, 8, 9, 11, 12, 13], 'ELG': [1, 11, 12, 13], 'QSO': [1, 8, 9, 11, 12, 13], 'BGS_ANY': [1, 13], 'BGS_BRIGHT': [1, 13]}
maskbits_dict = {'LRG': [], 'ELG': [1, 11, 12, 13], 'QSO': [1, 8, 9, 11, 12, 13], 'BGS_ANY': [1, 13], 'BGS_BRIGHT': [1, 13]}

randoms_counts_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_sys/randoms_stats/0.49.0/resolve/counts'
randoms_systematics_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_sys/randoms_stats/0.49.0/resolve/systematics'
target_densities_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_sys/density_maps/{}/resolve'.format(target_ver_str)

top_plot_dir = '/Users/rongpu/Documents/Work/DESI/imaging_sys/density_maps/{}/resolve'.format(target_ver_str)

dpi_dict = {64: 200, 128: 200, 256: 600, 512: 1600}
xsize_dict = {64: 8000, 128: 8000, 256: 12000, 512: 16000}
vrange_dict = {'BGS_ANY': {64: [800, 2000], 128: [650, 2150], 256: [200, 2600], 512: [-200, 3000]},
               'BGS_BRIGHT': {64: [500, 1200], 128: [350, 1350], 256: [200, 1500], 512: [-200, 1800]},
               'LRG': {64: [300, 900], 128: [200, 1000], 256: [100, 1100], 512: [-200, 1400]},
               'ELG': {64: [1200, 3600], 128: [1200, 3600], 256: [1100, 3700], 512: [600, 4200]},
               'QSO': {64: [150, 450], 128: [150, 450], 256: [0, 600], 512: [-200, 800]},
               }
# vrange_dict = {64: [0, 1200], 128: [-200, 1400], 256: [-600, 1800]}

min_pix_frac = 0.2  # minimum fraction of pixel area to be used

# for index in range(len(xnames)):
#     print(xnames[index], xlabels[index])


for target_class in ['BGS_ANY', 'BGS_BRIGHT', 'LRG', 'ELG', 'QSO']:

    if target_class=='LRG':
        lrgmask_str = '_lrgmask_v1'
    else:
        lrgmask_str = ''

    print(target_class)
    target_class = target_class.lower()

    maskbits = maskbits_dict[target_class.upper()]

    for nside in [64, 128, 256, 512]:

        npix = hp.nside2npix(nside)
        pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
        print(nside, 'Healpix size = {:.5f} sq deg'.format(pix_area))

        field = 'combined'

        plot_dir = os.path.join(top_plot_dir, '{}_{}_minobs_{}_maskbits_{}'.format(target_class, field, min_nobs, ''.join([str(tmp) for tmp in maskbits])+lrgmask_str))
        if not os.path.isdir(plot_dir):
            os.makedirs(plot_dir)

        for field in ['north', 'south']:

            density = Table.read(os.path.join(target_densities_dir, 'density_map_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(target_class, field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])+lrgmask_str)))
            maps = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])+lrgmask_str)))
            maps = maps[maps['n_randoms']>0]
            maps1 = Table.read(os.path.join(randoms_systematics_dir, 'systematics_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])+lrgmask_str)))
            maps1.remove_columns(['RA', 'DEC'])
            maps = join(maps, maps1, join_type='inner', keys='HPXPIXEL')
            maps = join(maps, density[['HPXPIXEL', 'n_targets']], join_type='outer', keys='HPXPIXEL').filled(0)

            # # Load stellar density map
            # stardens = np.load('/Users/rongpu/Documents/Data/desi_lrg_selection/dr7/healpix_maps/pixweight-dr7.1-0.22.0_stardens_{}_ring.npy'.format(nside))
            # maps['stardens'] = stardens[maps['HPXPIXEL']]
            # maps['stardens_log'] = np.log10(maps['stardens'])

            if field=='north':
                maps_north = maps.copy()
            else:
                maps_south = maps.copy()

        ########## Combine the two maps; proper handling of overlapping pixels ##########

        pix_overlap = np.intersect1d(maps_north['HPXPIXEL'], maps_south['HPXPIXEL'])
        mask = np.in1d(maps_north['HPXPIXEL'], pix_overlap)
        maps_overlap_north = maps_north[mask]
        maps_north = maps_north[~mask]
        mask = np.in1d(maps_south['HPXPIXEL'], pix_overlap)
        maps_overlap_south = maps_south[mask]
        maps_south = maps_south[~mask]

        maps_overlap_north.sort('HPXPIXEL')
        maps_overlap_south.sort('HPXPIXEL')

        maps_overlap = maps_overlap_south.copy()
        maps_overlap['n_targets'] = maps_overlap_north['n_targets'] + maps_overlap_south['n_targets']
        maps_overlap['FRACAREA'] = maps_overlap_north['FRACAREA'] + maps_overlap_south['FRACAREA']

        maps = vstack([maps_north, maps_south, maps_overlap])

        ######################################################################

        print(len(maps))

        area = np.sum(maps['FRACAREA'])*pix_area
        print('Area = {:.1f} sq deg'.format(area))

        mask = maps['FRACAREA']>min_pix_frac
        maps = maps[mask]

        maps['density'] = maps['n_targets'] / (pix_area * maps['FRACAREA'])

        plot_path = os.path.join(plot_dir, 'density_{}_{}{}.png'.format(target_class, nside, lrgmask_str))

        # if os.path.isfile(plot_path):
        #     continue

        map_values = np.zeros(npix)
        hp_mask = np.zeros(npix, dtype=bool)
        map_values[maps['HPXPIXEL']] = maps['density']
        hp_mask[maps['HPXPIXEL']] = True
        mplot = hp.ma(map_values)
        mplot.mask = ~hp_mask

        plt.figure(figsize=(9.7, 6))
        hp.mollview(mplot, min=vrange_dict[target_class.upper()][nside][0], max=vrange_dict[target_class.upper()][nside][1],
                    rot=(120, 0, 0), fig=1, xsize=xsize_dict[nside], title='{} NSIDE={}'.format(target_class.upper(), nside))
        plt.savefig(plot_path, dpi=dpi_dict[nside])
        plt.close()
