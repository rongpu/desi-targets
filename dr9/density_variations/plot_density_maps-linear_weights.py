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
maskbits_dict = {'LRG': [1, 8, 9, 11, 12, 13], 'ELG': [1, 11, 12, 13], 'QSO': [1, 8, 9, 11, 12, 13], 'BGS_ANY': [1, 13], 'BGS_BRIGHT': [1, 13]}

randoms_counts_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_systematics/randoms_stats/0.49.0/resolve/counts'
randoms_systematics_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_systematics/randoms_stats/0.49.0/resolve/systematics'
target_densities_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_systematics/density_maps/{}/resolve'.format(target_ver_str)
weights_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_systematics/randoms_stats/0.49.0/resolve/systematics/linear_weights_v0.1'

top_plot_dir = '/Users/rongpu/Documents/Work/DESI/imaging_systematics/density_maps/{}/resolve/linear_weights_v0.1'.format(target_ver_str)

dpi_dict = {64: 200, 128: 200, 256: 600, 512: 1600}
xsize_dict = {64: 8000, 128: 8000, 256: 12000, 512: 16000}
vrange_dict = {'BGS_ANY': {64: [800, 2000], 128: [650, 2150], 256: [400, 2400], 512: [0, 2800]},
               'BGS_BRIGHT': {64: [500, 1200], 128: [350, 1350], 256: [200, 1500], 512: [-200, 1800]},
               'LRG': {64: [300, 900], 128: [200, 1000], 256: [100, 1100], 512: [-100, 1300]},
               'ELG': {64: [1200, 3600], 128: [1200, 3600], 256: [1200, 3600], 512: [800, 4000]},
               'QSO': {64: [150, 450], 128: [150, 450], 256: [100, 500], 512: [0, 600]},
               }
# vrange_dict = {64: [0, 1200], 128: [-200, 1400], 256: [-600, 1800]}

min_pix_frac = 0.2  # minimum fraction of pixel area to be used

# for index in range(len(xnames)):
#     print(xnames[index], xlabels[index])


for target_class in ['BGS_ANY', 'BGS_BRIGHT', 'LRG', 'ELG', 'QSO']:

    print(target_class)
    target_class = target_class.lower()

    maskbits = maskbits_dict[target_class.upper()]

    for nside in [64, 128]:

        npix = hp.nside2npix(nside)
        pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
        print(nside, 'Healpix size = {:.5f} sq deg'.format(pix_area))

        field = 'combined'

        plot_dir = os.path.join(top_plot_dir, '{}_{}_minobs_{}_maskbits_{}'.format(target_class, field, min_nobs, ''.join([str(tmp) for tmp in maskbits])))
        if not os.path.isdir(plot_dir):
            os.makedirs(plot_dir)

        density_north = Table.read(os.path.join(target_densities_dir, 'density_map_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(target_class, 'north', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        density_south = Table.read(os.path.join(target_densities_dir, 'density_map_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(target_class, 'south', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        mask = (density_north['DEC']>32.375)
        density_north = density_north[mask]
        mask = ~np.in1d(density_south['HPXPIXEL'], density_north['HPXPIXEL'])
        density = vstack([density_north, density_south[mask]])

        weights_north = Table.read(os.path.join(weights_dir, 'linear_weights_v0.1_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(target_class, 'north', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        weights_south = Table.read(os.path.join(weights_dir, 'linear_weights_v0.1_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(target_class, 'south', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        mask = (weights_north['DEC']>32.375)
        weights_north = weights_north[mask]
        mask = ~np.in1d(weights_south['HPXPIXEL'], weights_north['HPXPIXEL'])
        weights = vstack([weights_north, weights_south[mask]])

        maps = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('north', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        maps = maps[maps['n_randoms']>0]
        maps_north = maps.copy()

        maps = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('south', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        maps = maps[maps['n_randoms']>0]
        maps_south = maps.copy()

        mask = (maps_north['DEC']>32.375)
        maps_north = maps_north[mask]
        mask = ~np.in1d(maps_south['HPXPIXEL'], maps_north['HPXPIXEL'])
        maps = vstack([maps_north, maps_south[mask]])

        maps = join(maps, density[['HPXPIXEL', 'n_targets']], join_type='outer', keys='HPXPIXEL').filled(0)
        maps = join(maps, weights[['HPXPIXEL', 'density_predict']], join_type='inner', keys='HPXPIXEL').filled(0)

        print(len(maps))

        area = np.sum(maps['FRACAREA'])*pix_area
        print('Area = {:.1f} sq deg'.format(area))

        mask = maps['FRACAREA']>min_pix_frac
        maps = maps[mask]

        maps['density'] = maps['n_targets'] / (pix_area * maps['FRACAREA'])
        maps['density_norm'] = maps['density'] / maps['density_predict']

        median_density = np.median(maps['density_predict'])

        ################### Plot weighted (and normalized) density map ###################

        plot_path = os.path.join(plot_dir, 'density_{}_{}_weighted_v0.1.png'.format(target_class, nside))

        if os.path.isfile(plot_path):
            continue

        map_values = np.zeros(npix)
        hp_mask = np.zeros(npix, dtype=bool)
        map_values[maps['HPXPIXEL']] = maps['density_norm'].copy()
        hp_mask[maps['HPXPIXEL']] = True
        mplot = hp.ma(map_values)
        mplot.mask = ~hp_mask

        plt.figure(figsize=(9.7, 6))
        hp.mollview(mplot, min=vrange_dict[target_class.upper()][nside][0]/median_density, max=vrange_dict[target_class.upper()][nside][1]/median_density,
                    rot=(120, 0, 0), fig=1, xsize=xsize_dict[nside], title='{} NSIDE={}'.format(target_class.upper(), nside))
        plt.savefig(plot_path, dpi=dpi_dict[nside])
        plt.close()

        ################### Plot "predicted" density map, i.e. weight map ###################

        plot_path = os.path.join(plot_dir, 'weight_map_{}_{}_v0.1.png'.format(target_class, nside))

        if os.path.isfile(plot_path):
            continue

        map_values = np.zeros(npix)
        hp_mask = np.zeros(npix, dtype=bool)
        map_values[maps['HPXPIXEL']] = maps['density_predict'].copy()
        hp_mask[maps['HPXPIXEL']] = True
        mplot = hp.ma(map_values)
        mplot.mask = ~hp_mask

        plt.figure(figsize=(9.7, 6))
        hp.mollview(mplot, min=vrange_dict[target_class.upper()][nside][0], max=vrange_dict[target_class.upper()][nside][1],
                    rot=(120, 0, 0), fig=1, xsize=xsize_dict[nside], title='{} NSIDE={}'.format(target_class.upper(), nside))
        plt.savefig(plot_path, dpi=dpi_dict[nside])
        plt.close()
