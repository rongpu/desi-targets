# Compute per-object linear weights for the LRG subsamples

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack
import fitsio

import yaml


def get_weights(cat, weights_path, bin_index):

    with open(weights_path, "r") as f:
        linear_coeffs = yaml.safe_load(f)

    weights = np.zeros(len(cat))

    regions = ['north', 'south']

    for region in regions:

        if region=='north' or region=='BASS+MzLS':
            photsys = 'N'
        else:
            photsys = 'S'
        
        bin_str = '{}_bin_{}'.format(region, bin_index)

        xnames_fit = list(linear_coeffs[bin_str].keys())
        xnames_fit.remove('intercept')

        # Assign nan weights to objects with invalid imaging properties
        # (their fraction should be negligibly small)
        if region==regions[0]:
            mask_bad = np.full(len(cat), False)
            for col in xnames_fit:
                mask_bad |= ~np.isfinite(cat[col])
            if np.sum(mask_bad)!=0:
                print('{} invalid objects'.format(np.sum(mask_bad)))
            weights[mask_bad] = np.nan

        # create array of coefficients, with the first coefficient being the intercept
        coeffs = np.array([linear_coeffs[bin_str]['intercept']]+[linear_coeffs[bin_str][xname] for xname in xnames_fit])

        mask = cat['PHOTSYS']==photsys
        mask &= (~mask_bad)

        data = np.column_stack([cat[mask][xname] for xname in xnames_fit])
        # create 2-D array of imaging properties, with the first columns being unity
        data1 = np.insert(data, 0, 1., axis=1)
        # wt = coeff0 + coeff1 * rand['EBV'] + coeff2 * rand['PSFSIZE_G'] + ...
        weights[mask] = np.dot(coeffs, data1.T)

    return weights


# # ##################### Example #####################

# # Prepare the catalog

# min_nobs = 2
# maskbits = sorted([1, 8, 9, 11, 12, 13])

# randoms_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'PHOTSYS',
#                    'GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z', 'PSFDEPTH_W1', 'PSFDEPTH_W2',
#                    'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 'EBV']

# cat = Table(fitsio.read('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-0.fits', columns=randoms_columns))
# # cat = Table(fitsio.read('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-0.fits', columns=randoms_columns,
# #                            rows=np.arange(int(1e6))))

# mask = (cat['NOBS_G']>=min_nobs) & (cat['NOBS_R']>=min_nobs) & (cat['NOBS_Z']>=min_nobs)
# cat = cat[mask]

# mask_clean = np.ones(len(cat), dtype=bool)
# for bit in maskbits:
#     mask_clean &= (cat['MASKBITS'] & 2**bit)==0
# cat = cat[mask_clean]

# with warnings.catch_warnings():
#     warnings.simplefilter("ignore")
#     cat['galdepth_gmag_ebv'] = -2.5*(np.log10((5/np.sqrt(cat['GALDEPTH_G'])))-9) - 3.214*cat['EBV']
#     cat['galdepth_rmag_ebv'] = -2.5*(np.log10((5/np.sqrt(cat['GALDEPTH_R'])))-9) - 2.165*cat['EBV']
#     cat['galdepth_zmag_ebv'] = -2.5*(np.log10((5/np.sqrt(cat['GALDEPTH_Z'])))-9) - 1.211*cat['EBV']
#     cat['psfdepth_w1mag_ebv'] = -2.5*(np.log10((5/np.sqrt(cat['PSFDEPTH_W1'])))-9) - 0.184*cat['EBV']
#     cat['psfdepth_w2mag_ebv'] = -2.5*(np.log10((5/np.sqrt(cat['PSFDEPTH_W2'])))-9) - 0.113*cat['EBV']

# # Get weights

# weights_path = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/imaging_weights/v1.1/main_lrg_linear_coeffs_pz.yaml'
# bin_index = 4

# weights = get_weights(cat, weights_path, bin_index)

