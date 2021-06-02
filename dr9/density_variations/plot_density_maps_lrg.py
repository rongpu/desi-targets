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

min_nobs = 2
maskbits = sorted([1, 8, 9, 11, 12, 13])

randoms_counts_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/randoms_stats/0.49.0/counts'
randoms_systematics_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/randoms_stats/0.49.0/systematics'
target_densities_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/unofficial/density_maps'

top_plot_dir = '/Users/rongpu/Documents/Work/DESI/imaging_systematics/density_maps/unofficial/noresolve'

target_name = 'sv3_lrg_all'
# target_name = 'sv3_lrg_lowdens'

dpi_dict = {64: 200, 128: 200, 256: 600}
vrange_dict = {64: [500, 1100], 128: [400, 1200], 256: [200, 1400]}
# vrange_dict = {64: [0, 1200], 128: [-200, 1400], 256: [-600, 1800]}

min_pix_frac = 0.2  # minimum fraction of pixel area to be used

# for index in range(len(xnames)):
#     print(xnames[index], xlabels[index])

# nside = 64
for nside in [64, 128, 256]:

    npix = hp.nside2npix(nside)
    pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
    print(nside, 'Healpix size = {:.5f} sq deg'.format(pix_area))

    for field in ['north', 'south', 'combined']:

        plot_dir = os.path.join(top_plot_dir, '{}_{}_minobs_{}_maskbits_{}'.format(target_name, field, min_nobs, ''.join([str(tmp) for tmp in maskbits])))
        if not os.path.isdir(plot_dir):
            os.makedirs(plot_dir)

        if field=='combined':

            density_north = Table.read(os.path.join(target_densities_dir, 'density_map_sv3_lrg_all_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('north', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
            density_south = Table.read(os.path.join(target_densities_dir, 'density_map_sv3_lrg_all_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('south', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
            # density_north = Table.read(os.path.join(target_densities_dir, 'density_map_sv3_lrg_lowdens_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('north', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
            # density_south = Table.read(os.path.join(target_densities_dir, 'density_map_sv3_lrg_lowdens_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('south', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
            mask = (density_north['dec']>32.375)
            density_north = density_north[mask]
            mask = ~np.in1d(density_south['hp_idx'], density_north['hp_idx'])
            density = vstack([density_north, density_south[mask]])
            density.rename_column('count', 'target_count')

            maps = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('north', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
            maps = maps[maps['count']>0]
            # maps1 = Table.read(os.path.join(randoms_systematics_dir, 'systematics_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('north', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
            # maps1.remove_columns(['ra', 'dec'])
            # maps = join(maps, maps1, join_type='inner', keys='hp_idx')
            maps_north = maps.copy()

            maps = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('south', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
            maps = maps[maps['count']>0]
            # maps1 = Table.read(os.path.join(randoms_systematics_dir, 'systematics_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('south', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
            # maps1.remove_columns(['ra', 'dec'])
            # maps = join(maps, maps1, join_type='inner', keys='hp_idx')
            maps_south = maps.copy()

            mask = (maps_north['dec']>32.375)
            maps_north = maps_north[mask]
            mask = ~np.in1d(maps_south['hp_idx'], maps_north['hp_idx'])
            maps = vstack([maps_north, maps_south[mask]])

        else:

            density = Table.read(os.path.join(target_densities_dir, 'density_map_sv3_lrg_all_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
            # density = Table.read(os.path.join(target_densities_dir, 'density_map_sv3_lrg_lowdens_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
            density.rename_column('count', 'target_count')

            maps = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
            maps = maps[maps['count']>0]
            # maps1 = Table.read(os.path.join(randoms_systematics_dir, 'systematics_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
            # maps1.remove_columns(['ra', 'dec'])
            # maps = join(maps, maps1, join_type='inner', keys='hp_idx')

        maps = join(maps, density[['hp_idx', 'NOBS_W1', 'target_count']], join_type='outer', keys='hp_idx').filled(0)

        print(len(maps))

        area = np.sum(maps['pix_frac'])*pix_area
        print('Area = {:.1f} sq deg'.format(area))

        mask = maps['pix_frac']>min_pix_frac
        maps = maps[mask]

        maps['density'] = maps['target_count'] / (pix_area * maps['pix_frac'])

        plot_path = os.path.join(plot_dir, 'density_{}_{}.png'.format(target_name, nside))

        # if os.path.isfile(plot_path):
        #     continue

        map_values = np.zeros(npix)
        hp_mask = np.zeros(npix, dtype=bool)
        map_values[maps['hp_idx']] = maps['density']
        hp_mask[maps['hp_idx']] = True
        mplot = hp.ma(map_values)
        mplot.mask = ~hp_mask

        plt.figure(figsize=(9.7, 6))
        hp.mollview(mplot, min=vrange_dict[nside][0], max=vrange_dict[nside][1],
                    rot=(120, 0, 0), fig=1, xsize=8000, title='{} NSIDE={}'.format(target_name, nside))
        plt.savefig(plot_path, dpi=dpi_dict[nside])
        plt.close()

