# NOT Separating DES and DECaLS footprints

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp
# from scipy import stats
# from IPython.display import Image

from sklearn.linear_model import LinearRegression

import yaml


randoms_counts_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/resolve/counts'
randoms_systematics_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/resolve/systematics'
target_densities_dir = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/density_maps/extended_lrg'
output_dir = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/imaging_weights/extended_lrg'

min_nobs = 2
maskbits = []
custom_mask_name = 'lrgmask_v1.1'
min_pix_frac = 0.6  # minimum fraction of pixel area to be used
nside = 256
include_ebv = True

mask_str = ''.join([str(tmp) for tmp in maskbits])
if custom_mask_name!='':
    mask_str += '_' + custom_mask_name

if include_ebv:
    xnames_fit = ['EBV', 'galdepth_gmag_ebv', 'galdepth_rmag_ebv', 'galdepth_zmag_ebv', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z']
    output_fn = 'extended_lrg_linear_coeffs_pz.yaml'
else:
    xnames_fit = ['galdepth_gmag_ebv', 'galdepth_rmag_ebv', 'galdepth_zmag_ebv', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z']
    output_fn = 'extended_lrg_linear_coeffs_pz_no_ebv.yaml'

coeff_dict = {}

for field in ['north', 'south']:

    for bin_index in range(1, 5):  # 4 bins

        print('bin {}'.format(bin_index))

        print('NSIDE = {}'.format(nside))

        npix = hp.nside2npix(nside)
        pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
        print('Healpix size = {:.5f} sq deg'.format(pix_area))

        density = Table(fitsio.read(os.path.join(target_densities_dir, 'density_map_extended_lrg_pz_bin_{}_{}_nside_{}_minobs_{}.fits'.format(bin_index, field, nside, min_nobs))))
        maps = Table(fitsio.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, mask_str))))
        maps = maps[maps['n_randoms']>0]
        maps1 = Table(fitsio.read(os.path.join(randoms_systematics_dir, 'systematics_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, mask_str))))
        maps1.remove_columns(['RA', 'DEC'])
        maps = join(maps, maps1, join_type='inner', keys='HPXPIXEL')
        maps = join(maps, density[['HPXPIXEL', 'n_targets']], join_type='outer', keys='HPXPIXEL').filled(0)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            maps['galdepth_gmag'] = -2.5*(np.log10((5/np.sqrt(maps['GALDEPTH_G'])))-9)
            maps['galdepth_rmag'] = -2.5*(np.log10((5/np.sqrt(maps['GALDEPTH_R'])))-9)
            maps['galdepth_zmag'] = -2.5*(np.log10((5/np.sqrt(maps['GALDEPTH_Z'])))-9)
            # maps['psfdepth_w1mag'] = -2.5*(np.log10((5/np.sqrt(maps['PSFDEPTH_W1'])))-9)
            # maps['psfdepth_w2mag'] = -2.5*(np.log10((5/np.sqrt(maps['PSFDEPTH_W2'])))-9)
            maps['galdepth_gmag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps['GALDEPTH_G'])))-9) - 3.214*maps['EBV']
            maps['galdepth_rmag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps['GALDEPTH_R'])))-9) - 2.165*maps['EBV']
            maps['galdepth_zmag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps['GALDEPTH_Z'])))-9) - 1.211*maps['EBV']
            # maps['psfdepth_w1mag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps['PSFDEPTH_W1'])))-9) - 0.184*maps['EBV']
            # maps['psfdepth_w2mag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps['PSFDEPTH_W2'])))-9) - 0.113*maps['EBV']

        mask = maps['FRACAREA']>min_pix_frac
        maps = maps[mask]

        # mask = maps['DEC']>-29  # Remove the southern part of DES
        # maps = maps[mask]

        # Remove pixels near the LMC
        ramin, ramax, decmin, decmax = 58, 110, -90, -56
        mask_remove = (maps['RA']>ramin) & (maps['RA']<ramax) & (maps['DEC']>decmin) & (maps['DEC']<decmax)
        maps = maps[~mask_remove]

        maps['density'] = maps['n_targets'] / (pix_area * maps['FRACAREA'])

        # Load stellar density map
        stardens = np.load('/global/cfs/cdirs/desi/users/rongpu/useful/healpix_maps/pixweight-dr7.1-0.22.0_stardens_{}_ring.npy'.format(nside))
        maps['stardens'] = stardens[maps['HPXPIXEL']]
        maps['stardens_log'] = np.log10(maps['stardens'])

        data = np.column_stack([maps[xname] for xname in xnames_fit])
        reg = LinearRegression().fit(data, maps['density'], sample_weight=maps['FRACAREA'])

        bin_str = '{}_bin_{}'.format(field, bin_index)
        coeff_dict[bin_str] = {}
        coeff_dict[bin_str]['intercept'] = float(reg.intercept_)
        for index, xname in enumerate(xnames_fit):
            coeff_dict[bin_str][xname] = float(reg.coef_[index])


with open(os.path.join(output_dir, output_fn), "w") as f:
    yaml.dump(coeff_dict, f)

