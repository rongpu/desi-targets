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

target_ver_str = '0.57.0'

min_nobs = 1
maskbits_dict = {'LRG': [1, 8, 9, 11, 12, 13], 'ELG': [1, 11, 12, 13], 'QSO': [1, 8, 9, 11, 12, 13], 'BGS_ANY': [1, 13], 'BGS_BRIGHT': [1, 13]}

randoms_counts_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_systematics/randoms_stats/0.49.0/resolve/counts'
randoms_systematics_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_systematics/randoms_stats/0.49.0/resolve/systematics'
target_densities_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_systematics/density_maps/{}/resolve'.format(target_ver_str)

top_plot_dir = '/Users/rongpu/Documents/Work/DESI/imaging_systematics/density_maps/{}/resolve'.format(target_ver_str)

dpi_dict = {64: 200, 128: 200, 256: 600, 512: 1600}
xsize_dict = {64: 8000, 128: 8000, 256: 12000, 512: 16000}
vrange_dict = {'SV3_BGS_ANY': {64: [800, 2000], 128: [650, 2150], 256: [400, 2400], 512: [0, 2800]},
               'SV3_BGS_BRIGHT': {64: [500, 1200], 128: [350, 1350], 256: [200, 1500], 512: [-200, 1800]},
               'SV3_LRG': {64: [500, 1100], 128: [400, 1200], 256: [300, 1300], 512: [100, 1500]},
               'SV3_ELG': {64: [1200, 3600], 128: [1200, 3600], 256: [1200, 3600], 512: [800, 4000]},
               'SV3_QSO': {64: [150, 450], 128: [150, 450], 256: [100, 500], 512: [0, 600]},
               }
# vrange_dict = {64: [0, 1200], 128: [-200, 1400], 256: [-600, 1800]}

min_pix_frac = 0.2  # minimum fraction of pixel area to be used

# for index in range(len(xnames)):
#     print(xnames[index], xlabels[index])


for target_class in ['SV3_BGS_ANY', 'SV3_BGS_BRIGHT', 'SV3_LRG', 'SV3_ELG', 'SV3_QSO']:

    print(target_class)
    target_class = target_class.lower()

    maskbits = maskbits_dict[target_class.upper()]

    for nside in [64, 128, 256, 512]:

        npix = hp.nside2npix(nside)
        pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
        print(nside, 'Healpix size = {:.5f} sq deg'.format(pix_area))

        field = 'combined'

        plot_dir = os.path.join(top_plot_dir, '{}_{}_minobs_{}_maskbits_{}'.format(target_class, field, min_nobs, ''.join([str(tmp) for tmp in maskbits])))
        if not os.path.isdir(plot_dir):
            os.makedirs(plot_dir)

        density_north = Table.read(os.path.join(target_densities_dir, 'density_map_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(target_class, 'north', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        density_south = Table.read(os.path.join(target_densities_dir, 'density_map_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(target_class, 'south', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        mask = (density_north['dec']>32.375)
        density_north = density_north[mask]
        mask = ~np.in1d(density_south['HPXPIXEL'], density_north['HPXPIXEL'])
        density = vstack([density_north, density_south[mask]])

        maps = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('north', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        maps = maps[maps['n_randoms']>0]
        # maps1 = Table.read(os.path.join(randoms_systematics_dir, 'systematics_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('north', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        # maps1.remove_columns(['ra', 'dec'])
        # maps = join(maps, maps1, join_type='inner', keys='HPXPIXEL')
        maps_north = maps.copy()

        maps = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('south', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        maps = maps[maps['n_randoms']>0]
        # maps1 = Table.read(os.path.join(randoms_systematics_dir, 'systematics_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('south', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        # maps1.remove_columns(['ra', 'dec'])
        # maps = join(maps, maps1, join_type='inner', keys='HPXPIXEL')
        maps_south = maps.copy()

        mask = (maps_north['dec']>32.375)
        maps_north = maps_north[mask]
        mask = ~np.in1d(maps_south['HPXPIXEL'], maps_north['HPXPIXEL'])
        maps = vstack([maps_north, maps_south[mask]])

        maps = join(maps, density[['HPXPIXEL', 'n_targets']], join_type='outer', keys='HPXPIXEL').filled(0)

        print(len(maps))

        area = np.sum(maps['FRACAREA'])*pix_area
        print('Area = {:.1f} sq deg'.format(area))

        mask = maps['FRACAREA']>min_pix_frac
        maps = maps[mask]

        maps['density'] = maps['n_targets'] / (pix_area * maps['FRACAREA'])

        plot_path = os.path.join(plot_dir, 'density_{}_{}.png'.format(target_class, nside))

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
