# source /global/common/software/desi/desi_environment.sh 22.2

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

plot_dir = '/global/cfs/cdirs/desicollab/users/rongpu/imaging_sys/density_maps/alternative_selections/ELG_LOP'

vrange_dict = {'BGS_ANY': {64: [800, 2000], 128: [650, 2150], 256: [200, 2600], 512: [-200, 3000]},
               'BGS_BRIGHT': {64: [500, 1200], 128: [350, 1350], 256: [200, 1500], 512: [-200, 1800]},
               'LRG': {64: [300, 900], 128: [200, 1000], 256: [100, 1100], 512: [-200, 1400]},
               'ELG': {64: [1200, 3600], 128: [1200, 3600], 256: [1100, 3700], 512: [600, 4200]},
               'ELG_LOP': {64: [1000, 2900], 128: [1000, 2900], 256: [900, 3000], 512: [500, 3400]},
               'QSO': {64: [150, 450], 128: [150, 450], 256: [0, 600], 512: [-200, 800]},
               }
# vrange_dict = {64: [0, 1200], 128: [-200, 1400], 256: [-600, 1800]}

min_pix_frac = 0.2  # minimum fraction of pixel area to be used

for ebv_str in ['desi_ebv', 'sfd_ebv', 'missing_elgs']:

    for selection in ['original', 'gmag', 'brighter', 'gmag_brighter']:

        for nside in [128]:

            npix = hp.nside2npix(nside)
            pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
            print(nside, 'Healpix size = {:.5f} sq deg'.format(pix_area))

            if ebv_str=='missing_elgs':
                map_path = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/misc/desi_targets_with_desi_ebv/density_map_alternative_{}_{}_{}.fits'.format(ebv_str, selection, nside)
            else:
                map_path = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/misc/desi_targets_with_desi_ebv/density_map_alternative_elgs_{}_{}_{}.fits'.format(ebv_str, selection, nside)

            plot_path = os.path.join(plot_dir, 'ELG_LOP_{}_{}_{}.png'.format(ebv_str, selection, nside))
            if os.path.isfile(plot_path):
                continue

            maps = Table(fitsio.read(map_path))
            print(len(maps))

            area = np.sum(maps['FRACAREA'])*pix_area
            print('Area = {:.1f} sq deg'.format(area))

            mask = maps['FRACAREA']>min_pix_frac
            maps = maps[mask]

            maps['density'] = maps['n_targets'] / (pix_area * maps['FRACAREA'])

            # Add the zero-density pixels for the "missing ELGs" map
            maps = maps['density', 'HPXPIXEL']
            if ebv_str=='missing_elgs':
                maps0 = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/misc/desi_targets_with_desi_ebv/density_map_alternative_elgs_{}_{}_{}.fits'.format('sfd_ebv', 'original', nside)))
                mask = maps0['FRACAREA']>min_pix_frac
                maps0 = maps0[mask]
                maps0 = maps0[['HPXPIXEL']]
            maps = join(maps, maps0, keys='HPXPIXEL', join_type='outer').filled(0.)

            if ebv_str in ['desi_ebv', 'sfd_ebv']:
                if selection in ['original', 'gmag']:
                    cmap, vmin, vmax = 'jet', vrange_dict['ELG_LOP'][nside][0], vrange_dict['ELG_LOP'][nside][1]
                elif selection=='gmag_brighter':
                    cmap, vmin, vmax = 'jet', 650 * 0.5, 650 * 1.5
                elif selection=='brighter':
                    cmap, vmin, vmax = 'jet', 570 * 0.5, 570 * 1.5
            elif ebv_str=='missing_elgs':
                if selection in ['original', 'gmag']:
                    cmap, vmin, vmax = 'bwr_r', -1000, 1000
                elif selection in ['brighter', 'gmag_brighter']:
                    cmap, vmin, vmax = 'bwr_r', -300, 300

            plot_map(nside, maps['density'], pix=maps['HPXPIXEL'], vmin=vmin, vmax=vmax,
                     title='ELG_LOP {} {} NSIDE={}'.format(ebv_str, selection, nside), save_path=plot_path, show=False, cmap=cmap, dpi=200, xsize=2000)

############################################# Difference between SFD-selected ELGs and DESI-selected ELGs #############################################

nside = 128
maps1 = Table(fitsio.read(f'/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/misc/desi_targets_with_desi_ebv/density_map_alternative_elgs_sfd_ebv_original_{nside}.fits'))
maps2 = Table(fitsio.read(f'/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/misc/desi_targets_with_desi_ebv/density_map_alternative_elgs_desi_ebv_original_{nside}.fits'))

pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
maps1['density'] = maps1['n_targets'] / (pix_area * maps1['FRACAREA'])
maps2['density'] = maps2['n_targets'] / (pix_area * maps2['FRACAREA'])

mask = maps1['FRACAREA']>0.2
print(np.sum(mask)/len(mask))
maps1 = maps1[mask]
mask = maps2['FRACAREA']>0.2
print(np.sum(mask)/len(mask))
maps2 = maps2[mask]

assert np.all(maps1['HPXPIXEL']==maps2['HPXPIXEL'])

plot_path = os.path.join(plot_dir, 'ELG_LOP_SFD_minus_ELG_LOP_DESI_{}.png'.format(nside))

plot_map(nside, maps1['density']-maps2['density'], pix=maps1['HPXPIXEL'], vmin=-1000, vmax=1000,
         title='ELG_LOP_SFD - ELG_LOP_DESI', save_path=plot_path, show=False, cmap='bwr')

# # Plot "Fraction of DESI footprint" vs density of missing ELGs
# density_list = np.arange(1000)
# fraction_list = np.zeros(len(density_list))
# for index, density in enumerate(density_list):
#     mask = maps1['density']-maps2['density']<-density
#     fraction_list[index] = np.sum(mask)/len(mask)
# plt.figure(figsize=(8, 6))
# plt.plot(density_list, fraction_list)
# plt.xlabel('density of "missing ELGs"(1/sq.deg.)')
# plt.ylabel('Fraction of DESI footprint')
# plt.grid(alpha=0.5)

