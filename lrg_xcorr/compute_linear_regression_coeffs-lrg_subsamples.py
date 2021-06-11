# NOT Separating DES and DECaLS footprints

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp
from scipy import stats
from IPython.display import Image

from sklearn.linear_model import LinearRegression

import yaml


randoms_counts_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_systematics/randoms_stats/0.49.0/resolve/counts'
randoms_systematics_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_systematics/randoms_stats/0.49.0/resolve/systematics'
target_densities_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/lrg_xcorr/density_maps/1.0.0/resolve'

min_nobs = 2
maskbits = [1, 8, 9, 11, 12, 13]

min_pix_frac = 0.6  # minimum fraction of pixel area to be used

nside = 256

output_dir = '/Users/rongpu/git/desi-targets/lrg_xcorr/data'

xnames_fit = ['EBV', 'psfdepth_w1mag_ebv', 'galdepth_gmag_ebv', 'galdepth_rmag_ebv', 'galdepth_zmag_ebv', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z']

####################################################### r-W1 bins #######################################################

coeff_dict = {}

for field in ['north', 'south']:

    for bin_index in range(1, 6):

        print('bin {}'.format(bin_index))

        print('NSIDE = {}'.format(nside))

        npix = hp.nside2npix(nside)
        pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
        print('Healpix size = {:.5f} sq deg'.format(pix_area))

        density = Table.read(os.path.join(target_densities_dir, 'density_map_lrg_rw1_bin_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(bin_index, field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        maps = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        maps = maps[maps['n_randoms']>0]
        maps1 = Table.read(os.path.join(randoms_systematics_dir, 'systematics_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        maps1.remove_columns(['RA', 'DEC'])
        maps = join(maps, maps1, join_type='inner', keys='HPXPIXEL')
        maps = join(maps, density[['HPXPIXEL', 'n_targets']], join_type='outer', keys='HPXPIXEL').filled(0)

        mask = maps['FRACAREA']>min_pix_frac
        maps = maps[mask]
        maps['density'] = maps['n_targets'] / (pix_area * maps['FRACAREA'])

        # Load stellar density map
        stardens = np.load('/Users/rongpu/Documents/Data/desi_lrg_selection/dr7/healpix_maps/pixweight-dr7.1-0.22.0_stardens_{}_ring.npy'.format(nside))
        maps['stardens'] = stardens[maps['HPXPIXEL']]
        maps['stardens_log'] = np.log10(maps['stardens'])

        data = np.column_stack([maps[xname] for xname in xnames_fit])
        reg = LinearRegression().fit(data, maps['density'], sample_weight=maps['FRACAREA'])

        bin_str = '{}_bin_{}'.format(field, bin_index)
        coeff_dict[bin_str] = {}
        coeff_dict[bin_str]['intercept'] = float(reg.intercept_)
        for index, xname in enumerate(xnames_fit):
            coeff_dict[bin_str][xname] = float(reg.coef_[index])

with open(os.path.join(output_dir, "main_lrg_linear_coeffs_rw1.yaml"), "w") as f:
    yaml.dump(coeff_dict, f)

####################################################### Photo-z bins #######################################################

coeff_dict = {}

for field in ['north', 'south']:

    for bin_index in range(1, 6):

        print('bin {}'.format(bin_index))

        print('NSIDE = {}'.format(nside))

        npix = hp.nside2npix(nside)
        pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
        print('Healpix size = {:.5f} sq deg'.format(pix_area))
    
        density = Table.read(os.path.join(target_densities_dir, 'density_map_lrg_pz_bin_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(bin_index, field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        maps = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        maps = maps[maps['n_randoms']>0]
        maps1 = Table.read(os.path.join(randoms_systematics_dir, 'systematics_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
        maps1.remove_columns(['RA', 'DEC'])
        maps = join(maps, maps1, join_type='inner', keys='HPXPIXEL')
        maps = join(maps, density[['HPXPIXEL', 'n_targets']], join_type='outer', keys='HPXPIXEL').filled(0)

        mask = maps['FRACAREA']>min_pix_frac
        maps = maps[mask]
        # mask = maps['DEC']>-30  # Remove the southern part of DES
        # maps = maps[mask]
        maps['density'] = maps['n_targets'] / (pix_area * maps['FRACAREA'])

        # Load stellar density map
        stardens = np.load('/Users/rongpu/Documents/Data/desi_lrg_selection/dr7/healpix_maps/pixweight-dr7.1-0.22.0_stardens_{}_ring.npy'.format(nside))
        maps['stardens'] = stardens[maps['HPXPIXEL']]
        maps['stardens_log'] = np.log10(maps['stardens'])

        data = np.column_stack([maps[xname] for xname in xnames_fit])
        reg = LinearRegression().fit(data, maps['density'], sample_weight=maps['FRACAREA'])

        bin_str = '{}_bin_{}'.format(field, bin_index)
        coeff_dict[bin_str] = {}
        coeff_dict[bin_str]['intercept'] = float(reg.intercept_)
        for index, xname in enumerate(xnames_fit):
            coeff_dict[bin_str][xname] = float(reg.coef_[index])

with open(os.path.join(output_dir, "main_lrg_linear_coeffs_pz.yaml"), "w") as f:
    yaml.dump(coeff_dict, f)

