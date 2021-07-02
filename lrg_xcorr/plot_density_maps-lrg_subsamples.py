# Subsample version 0.1

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
maskbits = [1, 8, 9, 11, 12, 13]

randoms_counts_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_systematics/randoms_stats/0.49.0/resolve/counts'
randoms_systematics_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_systematics/randoms_stats/0.49.0/resolve/systematics'
target_densities_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/lrg_xcorr/density_maps/1.0.0/resolve'

top_plot_dir = '/Users/rongpu/Documents/Work/DESI/lrg_xcorr/imaging_systematics/density_maps'

dpi_dict = {64: 200, 128: 200, 256: 600, 512: 1600}
xsize_dict = {64: 8000, 128: 8000, 256: 12000, 512: 16000}
vrange_dict = {64: [50, 170], 128: [30, 210], 256: [-30, 270]}
# vrange_dict = {64: [0, 1200], 128: [-200, 1400], 256: [-600, 1800]}

min_pix_frac = 0.2  # minimum fraction of pixel area to be used

# for index in range(len(xnames)):
#     print(xnames[index], xlabels[index])

############################## r-W1 bins ##############################

for bin_index in range(1, 6):

    # for nside in [64, 128, 256, 512]:
    for nside in [64, 128, 256]:

        npix = hp.nside2npix(nside)
        pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
        print(nside, 'Healpix size = {:.5f} sq deg'.format(pix_area))

        field = 'combined'

        plot_dir = os.path.join(top_plot_dir, 'lrg_rw1_{}_minobs_{}_maskbits_{}'.format(field, min_nobs, ''.join([str(tmp) for tmp in maskbits])))
        if not os.path.isdir(plot_dir):
            os.makedirs(plot_dir)

        plot_path = os.path.join(plot_dir, 'density_lrg_rw1_bin_{}_{}.png'.format(bin_index, nside))
        if os.path.isfile(plot_path):
            continue

        density_north = Table.read(os.path.join(target_densities_dir, 'density_map_lrg_rw1_bin_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(bin_index, 'north', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        density_south = Table.read(os.path.join(target_densities_dir, 'density_map_lrg_rw1_bin_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(bin_index, 'south', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        mask = (density_north['DEC']>32.375)
        density_north = density_north[mask]
        mask = ~np.in1d(density_south['HPXPIXEL'], density_north['HPXPIXEL'])
        density = vstack([density_north, density_south[mask]])

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

        print(len(maps))

        area = np.sum(maps['FRACAREA'])*pix_area
        print('Area = {:.1f} sq deg'.format(area))

        mask = maps['FRACAREA']>min_pix_frac
        maps = maps[mask]

        maps['density'] = maps['n_targets'] / (pix_area * maps['FRACAREA'])

        map_values = np.zeros(npix)
        hp_mask = np.zeros(npix, dtype=bool)
        map_values[maps['HPXPIXEL']] = maps['density']
        hp_mask[maps['HPXPIXEL']] = True
        mplot = hp.ma(map_values)
        mplot.mask = ~hp_mask

        plt.figure(figsize=(9.7, 6))
        hp.mollview(mplot, min=vrange_dict[nside][0], max=vrange_dict[nside][1],
                    rot=(120, 0, 0), fig=1, xsize=xsize_dict[nside], title='LRG r-W1 bin {} NSIDE={}'.format(bin_index, nside))
        plt.savefig(plot_path, dpi=dpi_dict[nside])
        plt.close()


############################## photo-z bins ##############################

for bin_index in range(1, 6):

    # for nside in [64, 128, 256, 512]:
    for nside in [64, 128, 256]:

        npix = hp.nside2npix(nside)
        pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
        print(nside, 'Healpix size = {:.5f} sq deg'.format(pix_area))

        field = 'combined'

        plot_dir = os.path.join(top_plot_dir, 'lrg_pz_{}_minobs_{}_maskbits_{}'.format(field, min_nobs, ''.join([str(tmp) for tmp in maskbits])))
        if not os.path.isdir(plot_dir):
            os.makedirs(plot_dir)

        plot_path = os.path.join(plot_dir, 'density_lrg_pz_bin_{}_{}.png'.format(bin_index, nside))
        if os.path.isfile(plot_path):
            continue

        density_north = Table.read(os.path.join(target_densities_dir, 'density_map_lrg_pz_bin_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(bin_index, 'north', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        density_south = Table.read(os.path.join(target_densities_dir, 'density_map_lrg_pz_bin_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(bin_index, 'south', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        mask = (density_north['DEC']>32.375)
        density_north = density_north[mask]
        mask = ~np.in1d(density_south['HPXPIXEL'], density_north['HPXPIXEL'])
        density = vstack([density_north, density_south[mask]])

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

        print(len(maps))

        area = np.sum(maps['FRACAREA'])*pix_area
        print('Area = {:.1f} sq deg'.format(area))

        mask = maps['FRACAREA']>min_pix_frac
        maps = maps[mask]

        maps['density'] = maps['n_targets'] / (pix_area * maps['FRACAREA'])

        map_values = np.zeros(npix)
        hp_mask = np.zeros(npix, dtype=bool)
        map_values[maps['HPXPIXEL']] = maps['density']
        hp_mask[maps['HPXPIXEL']] = True
        mplot = hp.ma(map_values)
        mplot.mask = ~hp_mask

        plt.figure(figsize=(9.7, 6))
        hp.mollview(mplot, min=vrange_dict[nside][0], max=vrange_dict[nside][1],
                    rot=(120, 0, 0), fig=1, xsize=xsize_dict[nside], title='LRG photo-z bin {} NSIDE={}'.format(bin_index, nside))
        plt.savefig(plot_path, dpi=dpi_dict[nside])
        plt.close()
