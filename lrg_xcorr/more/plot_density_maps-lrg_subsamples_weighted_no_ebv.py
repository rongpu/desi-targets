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

weighted = True

min_nobs = 2
maskbits = []
custom_mask_name = 'lrgmask_v1.1'

mask_str = ''.join([str(tmp) for tmp in maskbits])
if custom_mask_name!='':
    mask_str += '_' + custom_mask_name

randoms_counts_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/resolve/counts'
randoms_systematics_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/resolve/systematics'
target_densities_dir = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/density_maps/main_lrg'
top_plot_dir = '/global/cfs/cdirs/desi/users/rongpu/lrg_xcorr/imaging_sys/density_maps/main_lrg'

if weighted:
    target_densities_dir = os.path.join(target_densities_dir, 'linear_weights')
    top_plot_dir = os.path.join(top_plot_dir, 'linear_weights')

dpi_dict = {64: 200, 128: 200, 256: 600, 512: 1600}
xsize_dict = {64: 8000, 128: 8000, 256: 12000, 512: 16000}
vrange_dict = {64: 8, 128: 12, 256: 18}
vcenter_dict = {1: 81.0, 2: 147.5, 3: 164.2, 4: 149.4}  # average density in DECam imaging
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

        plot_dir = os.path.join(top_plot_dir, 'lrg_pz_{}_minobs_{}_maskbits_{}'.format(field, min_nobs, mask_str))
        if not os.path.isdir(plot_dir):
            os.makedirs(plot_dir)

        plot_path = os.path.join(plot_dir, 'density_lrg_pz_bin_{}_{}.png'.format(bin_index, nside))
        if weighted:
            plot_path = plot_path.replace('.png', '-lw_no_ebv.png')
        if os.path.isfile(plot_path):
            continue

        if weighted:
            weighted_str = '-lw_no_ebv'
        else:
            weighted_str = ''

        for field in ['north', 'south']:

            density = Table(fitsio.read(os.path.join(target_densities_dir, 'density_map_lrg_pz_bin_{}_{}_nside_{}_minobs_{}{}.fits'.format(bin_index, field, nside, min_nobs, weighted_str))))
            maps = Table(fitsio.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, mask_str))))
            maps = maps[maps['n_randoms']>0]
            maps1 = Table(fitsio.read(os.path.join(randoms_systematics_dir, 'systematics_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, mask_str))))
            maps1.remove_columns(['RA', 'DEC'])
            maps = join(maps, maps1, join_type='inner', keys='HPXPIXEL')
            maps = join(maps, density[['HPXPIXEL', 'n_targets']], join_type='outer', keys='HPXPIXEL').filled(0)

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

        area = np.sum(maps['FRACAREA'])*pix_area
        print('Area = {:.1f} sq deg'.format(area))

        mask = maps['FRACAREA']>min_pix_frac
        maps = maps[mask]

        maps['density'] = maps['n_targets'] / (pix_area * maps['FRACAREA'])

        vcenter = vcenter_dict[bin_index]
        vrange = int(np.sqrt(vcenter)*vrange_dict[nside])
        vmin, vmax = vcenter - vrange, vcenter + vrange
        if weighted:
            vmin, vmax = vmin/vcenter, vmax/vcenter

        plot_map(nside, maps['HPXPIXEL'], maps['density'],
                 vmin=vmin, vmax=vmax, xsize=xsize_dict[nside], dpi=dpi_dict[nside],
                 title='LRG photo-z bin {} NSIDE={}'.format(bin_index, nside), save_path=plot_path, show=False)
