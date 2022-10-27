# version 1.1
# This script assumes that the "resolve" catalogs are used
# The overlaping pixels are combined by weighted averaging

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

weighted = False

min_nobs = 2
# maskbits = [1, 8, 9, 11, 12, 13]

maskbits = []
apply_lrgmask = True
if apply_lrgmask:
    lrgmask_str = '_lrgmask_v1'
else:
    lrgmask_str = ''

randoms_counts_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_sys/randoms_stats/0.49.0/resolve/counts'
randoms_systematics_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_sys/randoms_stats/0.49.0/resolve/systematics'
target_densities_dir = '/Users/rongpu/Documents/Data/lrg_xcorr/density_maps/1.0.0/resolve/v1.1'

top_plot_dir = '/Users/rongpu/Documents/Work/DESI/lrg_xcorr/imaging_sys/density_maps/v1.1'

if weighted:
    target_densities_dir = os.path.join(target_densities_dir, 'linear_weights')
    top_plot_dir = os.path.join(top_plot_dir, 'linear_weights')

dpi_dict = {64: 200, 128: 200, 256: 600, 512: 1600}
xsize_dict = {64: 8000, 128: 8000, 256: 12000, 512: 16000}
vrange_dict = {64: 6, 128: 10, 256: 15}
vcenter_dict = {1: 85, 2: 150, 3: 165, 4: 150}
# vrange_dict = {64: [0, 1200], 128: [-200, 1400], 256: [-600, 1800]}

min_pix_frac = 0.2  # minimum fraction of pixel area to be used

# for index in range(len(xnames)):
#     print(xnames[index], xlabels[index])

for bin_index in range(1, 5):

    for nside in [64, 128, 256]:

        npix = hp.nside2npix(nside)
        pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
        print(nside, 'Healpix size = {:.5f} sq deg'.format(pix_area))

        field = 'combined'

        plot_dir = os.path.join(top_plot_dir, 'lrg_pz_{}_minobs_{}_maskbits_{}'.format(field, min_nobs, ''.join([str(tmp) for tmp in maskbits])))
        if not os.path.isdir(plot_dir):
            os.makedirs(plot_dir)

        plot_path = os.path.join(plot_dir, 'density_lrg_pz_bin_{}_{}.png'.format(bin_index, nside))
        if weighted:
            plot_path = plot_path.replace('.png', '-lw.png')
        if os.path.isfile(plot_path):
            continue

        if weighted:
            weighted_str = '-lw'
        else:
            weighted_str = ''
        # density_north = Table.read(os.path.join(target_densities_dir, 'density_map_lrg_pz_bin_{}_{}_nside_{}_minobs_{}_maskbits_{}{}.fits'.format(bin_index, 'north', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]), weighted_str)))
        # density_south = Table.read(os.path.join(target_densities_dir, 'density_map_lrg_pz_bin_{}_{}_nside_{}_minobs_{}_maskbits_{}{}.fits'.format(bin_index, 'south', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]), weighted_str)))
        density_north = Table.read(os.path.join(target_densities_dir, 'density_map_lrg_pz_bin_{}_{}_nside_{}_minobs_{}{}.fits'.format(bin_index, 'north', nside, min_nobs, weighted_str)))
        maps_north = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('north', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])+lrgmask_str)))
        maps_north = maps_north[maps_north['n_randoms']>0]
        maps_north = join(maps_north, density_north[['HPXPIXEL', 'n_targets']], join_type='outer', keys='HPXPIXEL').filled(0)

        density_south = Table.read(os.path.join(target_densities_dir, 'density_map_lrg_pz_bin_{}_{}_nside_{}_minobs_{}{}.fits'.format(bin_index, 'south', nside, min_nobs, weighted_str)))
        maps_south = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format('south', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])+lrgmask_str)))
        maps_south = maps_south[maps_south['n_randoms']>0]
        maps_south = join(maps_south, density_south[['HPXPIXEL', 'n_targets']], join_type='outer', keys='HPXPIXEL').filled(0)

        # Combine north and south
        pix_overlap = np.intersect1d(maps_north['HPXPIXEL'], maps_south['HPXPIXEL'])
        maps_north_overlap = maps_north[np.in1d(maps_north['HPXPIXEL'], pix_overlap)].copy()
        maps_north_overlap.sort('HPXPIXEL')
        maps_south_overlap = maps_south[np.in1d(maps_south['HPXPIXEL'], pix_overlap)].copy()
        maps_south_overlap.sort('HPXPIXEL')
        maps_overlap = maps_north_overlap.copy()
        maps_overlap['n_targets'] = maps_north_overlap['n_targets'] + maps_south_overlap['n_targets']
        maps_overlap['n_randoms'] = maps_north_overlap['n_randoms'] + maps_south_overlap['n_randoms']
        maps_overlap['FRACAREA'] = maps_north_overlap['FRACAREA'] + maps_south_overlap['FRACAREA']

        maps_north = maps_north[~np.in1d(maps_north['HPXPIXEL'], pix_overlap)]
        maps_south = maps_south[~np.in1d(maps_south['HPXPIXEL'], pix_overlap)]
        maps = vstack([maps_north, maps_south, maps_overlap], join_type='exact')

        mask = maps['FRACAREA']>min_pix_frac
        maps = maps[mask]

        maps['density'] = maps['n_targets'] / (pix_area * maps['FRACAREA'])

        map_values = np.zeros(npix)
        hp_mask = np.zeros(npix, dtype=bool)
        map_values[maps['HPXPIXEL']] = maps['density']
        hp_mask[maps['HPXPIXEL']] = True
        mplot = hp.ma(map_values)
        mplot.mask = ~hp_mask

        vcenter = vcenter_dict[bin_index]
        vrange = int(np.sqrt(vcenter)*vrange_dict[nside])
        vmin, vmax = vcenter - vrange, vcenter + vrange
        if weighted:
            vmin, vmax = vmin/vcenter, vmax/vcenter

        # vrange = vrange_dict[nside]
        plt.figure(figsize=(9.7, 6))
        hp.mollview(mplot, min=vmin, max=vmax,
                    rot=(120, 0, 0), fig=1, xsize=xsize_dict[nside], title='LRG photo-z bin {} NSIDE={}'.format(bin_index, nside))
        plt.savefig(plot_path, dpi=dpi_dict[nside])
        plt.close()
