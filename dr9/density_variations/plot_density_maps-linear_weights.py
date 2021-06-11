# NSIDE=512 for weights and downsize to NSIDE=64 for plotting

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import yaml

import healpy as hp
from sklearn.linear_model import LinearRegression

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

top_plot_dir = '/Users/rongpu/Documents/Work/DESI/imaging_systematics/density_maps/{}/resolve/linear_weights'.format(target_ver_str)
weights_dir ='/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_systematics/linear_weights'

dpi_dict = {64: 200, 128: 200, 256: 600, 512: 1600}
xsize_dict = {64: 8000, 128: 8000, 256: 12000, 512: 16000}
vrange_dict = {'BGS_ANY': {64: [800, 2000], 128: [650, 2150], 256: [400, 2400], 512: [0, 2800]},
               'BGS_BRIGHT': {64: [500, 1200], 128: [350, 1350], 256: [200, 1500], 512: [-200, 1800]},
               'LRG': {64: [300, 900], 128: [200, 1000], 256: [100, 1100], 512: [-100, 1300]},
               'ELG': {64: [1200, 3600], 128: [1200, 3600], 256: [1200, 3600], 512: [800, 4000]},
               'QSO': {64: [150, 450], 128: [150, 450], 256: [100, 500], 512: [0, 600]},
               }
# vrange_dict = {64: [0, 1200], 128: [-200, 1400], 256: [-600, 1800]}

min_pix_frac = 0.6  # minimum fraction of pixel area to be used

nside = 512
nside_plot = 64

# for index in range(len(xnames)):
#     print(xnames[index], xlabels[index])


for target_class in ['BGS_ANY', 'BGS_BRIGHT', 'LRG', 'ELG', 'QSO']:

    print(target_class)
    target_class = target_class.lower()

    maskbits = maskbits_dict[target_class.upper()]

    with open(os.path.join(weights_dir, "main_{}_linear_coeffs_v0.1.yaml".format(target_class)), "r") as f:
        linear_coeffs = yaml.safe_load(f)

    npix = hp.nside2npix(nside_weight)
    pix_area = hp.pixelfunc.nside2pixarea(nside_weight, degrees=True)
    print(nside_weight, 'Healpix size = {:.5f} sq deg'.format(pix_area))

    field = 'combined'

    plot_dir = os.path.join(top_plot_dir, '{}_{}_minobs_{}_maskbits_{}'.format(target_class, field, min_nobs, ''.join([str(tmp) for tmp in maskbits])))
    if not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)

    maps_dict = {}
    for field in ['north', 'south']:

        density = Table.read(os.path.join(target_densities_dir, 'density_map_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(target_class, field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        maps = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        maps = maps[maps['n_randoms']>0]
        maps1 = Table.read(os.path.join(randoms_systematics_dir, 'systematics_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        maps1.remove_columns(['RA', 'DEC'])
        maps = join(maps, maps1, join_type='inner', keys='HPXPIXEL')
        maps = join(maps, density[['HPXPIXEL', 'n_targets']], join_type='outer', keys='HPXPIXEL').filled(0)

        mask = maps['FRACAREA']>min_pix_frac
        maps = maps[mask]
        mask = maps['DEC']>-30  # Remove the southern part of DES
        maps = maps[mask]
        maps['density'] = maps['n_targets'] / (pix_area * maps['FRACAREA'])

        # Load stellar density map
        stardens = np.load('/Users/rongpu/Documents/Data/desi_lrg_selection/dr7/healpix_maps/pixweight-dr7.1-0.22.0_stardens_{}_ring.npy'.format(nside))
        maps['stardens'] = stardens[maps['HPXPIXEL']]
        maps['stardens_log'] = np.log10(maps['stardens'])

        maps_dict[field] = maps.copy()

    # Combine the maps
    mask = (maps_dict['north']['DEC']>32.375)
    mask1 = ~np.in1d(maps_dict['south']['HPXPIXEL'], maps_dict['north']['HPXPIXEL'][mask])
    maps = vstack([maps_dict['north'][mask], maps_dict['south'][mask1]]).copy()
    print(len(maps))

    area = np.sum(maps['FRACAREA'])*pix_area
    print('Area = {:.1f} sq deg'.format(area))

    mask = maps['FRACAREA']>min_pix_frac
    maps = maps[mask]

    maps['density'] = maps['n_targets'] / (pix_area * maps['FRACAREA'])

    ############################################ Apply weights #####################################################
    
    for field in ['north', 'south']:
        
        xnames_fit = list(linear_coeffs[field].keys())
        xnames_fit.remove('intercept')

        reg = LinearRegression()
        reg.intercept_ = linear_coeffs[field]['intercept']
        reg.coef_ = np.array([linear_coeffs[field][xname] for xname in xnames_fit])
        
        data = np.column_stack([maps[field][xname] for xname in xnames_fit])
        maps[field]['density_predict'] = reg.predict(data)

    ########################################################################################################

    plot_path = os.path.join(plot_dir, 'density_{}_{}.png'.format(target_class, nside_plot))

    # if os.path.isfile(plot_path):
    #     continue

    map_values = np.zeros(npix)
    hp_mask = np.zeros(npix, dtype=bool)
    map_values[maps['HPXPIXEL']] = maps['density']
    hp_mask[maps['HPXPIXEL']] = True
    mplot = hp.ma(map_values)
    mplot.mask = ~hp_mask

    plt.figure(figsize=(9.7, 6))
    hp.mollview(mplot, min=vrange_dict[target_class.upper()][nside_plot][0], max=vrange_dict[target_class.upper()][nside_plot][1],
                rot=(120, 0, 0), fig=1, xsize=xsize_dict[nside_plot], title='{} NSIDE={}'.format(target_class.upper(), nside_plot))
    plt.savefig(plot_path, dpi=dpi_dict[nside_plot])
    plt.close()
