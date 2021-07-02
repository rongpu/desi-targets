from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table
import fitsio

import yaml


def get_randoms_weights(randoms, weights_path, separate_des=False):

    with open(weights_path, "r") as f:
        linear_coeffs = yaml.safe_load(f)

    weights = np.zeros(len(randoms))

    if separate_des==False:  # Same linear coefficients for DECaLS and DES

        for field in ['north', 'south']:

            if field=='south':
                photsys = 'S'
            elif field=='north':
                photsys = 'N'

            xnames_fit = list(linear_coeffs[field].keys())
            xnames_fit.remove('intercept')

            # create array of coefficients, with the first coefficient being the intercept
            coeffs = np.array([linear_coeffs[field]['intercept']]+[linear_coeffs[field][xname] for xname in xnames_fit])

            mask = randoms['PHOTSYS']==photsys

            # Assign zero weights to randoms with invalid imaging properties
            # (their fraction should be negligibly small)
            mask_bad = np.full(len(randoms), False)
            for col in xnames_fit:
                mask_bad |= ~np.isfinite(randoms[col])
            if np.sum(mask_bad)!=0:
                print('{} invalid randoms'.format(np.sum(mask_bad)))
            mask &= (~mask_bad)

            data = np.column_stack([randoms[mask][xname] for xname in xnames_fit])
            # create 2-D array of imaging properties, with the first columns being unity
            data1 = np.insert(data, 0, 1., axis=1)
            # wt = coeff0 + coeff1 * rand['EBV'] + coeff2 * rand['PSFSIZE_G'] + ...
            weights[mask] = np.dot(coeffs, data1.T)

    else:  # Different linear coefficients for DECaLS and DES

        for region in ['BASS+MzLS', 'DECaLS', 'DES']:

            if region=='BASS+MzLS':
                field = 'north'
                photsys = 'N'
            else:
                field = 'south'
                photsys = 'S'

            xnames_fit = list(linear_coeffs[region].keys())
            xnames_fit.remove('intercept')

            # create array of coefficients, with the first coefficient being the intercept
            coeffs = np.array([linear_coeffs[region]['intercept']]+[linear_coeffs[region][xname] for xname in xnames_fit])

            mask = randoms['PHOTSYS']==photsys

            if field=='south':
                # See https://github.com/desihub/desitarget/blob/f39455769628b7982fc18c1e2668a6d1161a3e87/py/desitarget/cuts.py#L1874
                is_des = (randoms['NOBS_G'] > 4) & (randoms['NOBS_R'] > 4) & (randoms['NOBS_Z'] > 4) \
                    & ((randoms['RA'] >= 320) | (randoms['RA'] <= 100)) &  (randoms['DEC'] <= 10)
                if region=='DECaLS':
                    mask &= (~is_des)
                elif region=='DES':
                    mask &= is_des

            # Assign zero weights to randoms with invalid imaging properties
            # (their fraction should be negligibly small)
            mask_bad = np.full(len(randoms), False)
            for col in xnames_fit:
                mask_bad |= ~np.isfinite(randoms[col])
            if np.sum(mask_bad)!=0:
                print('{} invalid randoms'.format(np.sum(mask_bad)))
            mask &= (~mask_bad)

            data = np.column_stack([randoms[mask][xname] for xname in xnames_fit])
            # create 2-D array of imaging properties, with the first columns being unity
            data1 = np.insert(data, 0, 1., axis=1)
            # wt = coeff0 + coeff1 * rand['EBV'] + coeff2 * rand['PSFSIZE_G'] + ...
            weights[mask] = np.dot(coeffs, data1.T)

    return weights


# Prepare the randoms

min_nobs = 1  # Minimum NOBS_G,R,Z; default is 1 for DESI targets
maskbits = sorted([1, 8, 9, 11, 12, 13])

randoms_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'PHOTSYS',
                   'GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z', 
                   'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'PSFDEPTH_W1', 'PSFDEPTH_W2',
                   'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 'EBV']

# randoms = Table(fitsio.read('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-0.fits', columns=randoms_columns))
randoms = Table(fitsio.read('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-0.fits', columns=randoms_columns,
                           rows=np.arange(int(1e6))))

mask = (randoms['NOBS_G']>=min_nobs) & (randoms['NOBS_R']>=min_nobs) & (randoms['NOBS_Z']>=min_nobs)
randoms = randoms[mask]

mask_clean = np.ones(len(randoms), dtype=bool)
for bit in maskbits:
    mask_clean &= (randoms['MASKBITS'] & 2**bit)==0
randoms = randoms[mask_clean]

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    randoms['galdepth_gmag_ebv'] = -2.5*(np.log10((5/np.sqrt(randoms['GALDEPTH_G'])))-9) - 3.214*randoms['EBV']
    randoms['galdepth_rmag_ebv'] = -2.5*(np.log10((5/np.sqrt(randoms['GALDEPTH_R'])))-9) - 2.165*randoms['EBV']
    randoms['galdepth_zmag_ebv'] = -2.5*(np.log10((5/np.sqrt(randoms['GALDEPTH_Z'])))-9) - 1.211*randoms['EBV']
    randoms['psfdepth_gmag_ebv'] = -2.5*(np.log10((5/np.sqrt(randoms['PSFDEPTH_G'])))-9) - 3.214*randoms['EBV']
    randoms['psfdepth_rmag_ebv'] = -2.5*(np.log10((5/np.sqrt(randoms['PSFDEPTH_R'])))-9) - 2.165*randoms['EBV']
    randoms['psfdepth_zmag_ebv'] = -2.5*(np.log10((5/np.sqrt(randoms['PSFDEPTH_Z'])))-9) - 1.211*randoms['EBV']
    randoms['psfdepth_w1mag_ebv'] = -2.5*(np.log10((5/np.sqrt(randoms['PSFDEPTH_W1'])))-9) - 0.184*randoms['EBV']
    randoms['psfdepth_w2mag_ebv'] = -2.5*(np.log10((5/np.sqrt(randoms['PSFDEPTH_W2'])))-9) - 0.113*randoms['EBV']

# Get weights
# weights_path = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/linear_weights/main_v0.1/main_lrg_linear_coeffs_v0.1.yaml'
weights_path = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/linear_weights/main_v0.1/main_qso_linear_coeffs_separate_des_v0.1.yaml'

if '_separate_des_' in os.path.basename(weights_path):
    separate_des = True
else:
    separate_des = False
print("Separate DES weights:", separate_des)

weights = get_randoms_weights(randoms, weights_path, separate_des=separate_des)

