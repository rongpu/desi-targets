# Different linear models for DES and DECaLS

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

sys.path.append(os.path.expanduser('~/git/desi-targets/useful'))
from isdes import get_isdes


randoms_counts_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_systematics/randoms_stats/0.49.0/resolve/counts'
randoms_systematics_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_systematics/randoms_stats/0.49.0/resolve/systematics'
target_densities_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_systematics/density_maps/1.0.0/resolve'

min_nobs = 1
maskbits_dict = {'LRG': [1, 8, 9, 11, 12, 13], 'ELG': [1, 11, 12, 13], 'QSO': [1, 8, 9, 11, 12, 13], 'BGS_ANY': [1, 13], 'BGS_BRIGHT': [1, 13]}

min_pix_frac = 0.6  # minimum fraction of pixel area to be used

nside = 256

output_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_systematics/linear_weights'

xnames_fit_dict = {'BGS_ANY':['EBV', 'galdepth_rmag_ebv', 'PSFSIZE_R'],
                   'BGS_BRIGHT':['EBV', 'galdepth_rmag_ebv', 'PSFSIZE_R'],
                   'LRG':['EBV', 'psfdepth_w1mag_ebv', 'galdepth_gmag_ebv', 'galdepth_rmag_ebv', 'galdepth_zmag_ebv', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z'],
                   'ELG':['EBV', 'galdepth_gmag_ebv', 'galdepth_rmag_ebv', 'galdepth_zmag_ebv', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z'],
                   'QSO':['EBV', 'psfdepth_w1mag_ebv', 'psfdepth_w2mag_ebv', 'psfdepth_gmag_ebv', 'psfdepth_rmag_ebv', 'psfdepth_zmag_ebv', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z']
                   }

print('NSIDE = {}'.format(nside))

npix = hp.nside2npix(nside)
pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
print('Healpix size = {:.5f} sq deg'.format(pix_area))

for target_class in ['BGS_ANY', 'BGS_BRIGHT', 'LRG', 'ELG', 'QSO']:
    
    print(target_class)
    target_class = target_class.lower()

    maskbits = maskbits_dict[target_class.upper()]
    xnames_fit = xnames_fit_dict[target_class.upper()]

    coeff_dict = {}

    for region in ['BASS+MzLS', 'DECaLS', 'DES']:

        if region=='BASS+MzLS':
            field = 'north'
        else:
            field = 'south'

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

        if field=='south':
            isdes = get_isdes(maps['RA'], maps['DEC'], nside)
            if region=='DES':
                maps = maps[isdes]
            else:
                maps = maps[~isdes]

        data = np.column_stack([maps[xname] for xname in xnames_fit])
        reg = LinearRegression().fit(data, maps['density'], sample_weight=maps['FRACAREA'])

        coeff_dict[region] = {}
        coeff_dict[region]['intercept'] = float(reg.intercept_)
        for index, xname in enumerate(xnames_fit):
            coeff_dict[region][xname] = float(reg.coef_[index])

    with open(os.path.join(output_dir, "main_{}_linear_coeffs_separate_des_v0.1.yaml".format(target_class)), "w") as f:
        yaml.dump(coeff_dict, f)

