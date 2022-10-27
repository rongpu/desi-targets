# Subsample v1.0

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio


def apply_mask(cat, min_nobs, maskbits):

    mask = (cat['NOBS_G']>=min_nobs) & (cat['NOBS_R']>=min_nobs) & (cat['NOBS_Z']>=min_nobs)

    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    # print(np.sum(~mask_clean)/len(mask_clean))

    mask &= mask_clean

    return mask


min_nobs = 1
maskbits = []

basic_columns = ['TARGETID', 'RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'PHOTSYS']
photom_columns = ['EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z',
                  'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'FLUX_W1',
                  'MW_TRANSMISSION_W1', 'FIBERFLUX_Z']

columns_to_keep = ['TARGETID', 'RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'PHOTSYS', 'Z_PHOT_MEDIAN', 'pz_bin']

cat_stack = []
for field in ['north', 'south']:

    cat0 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve/dr9_lrg_{}_1.0.0_basic.fits'.format(field), columns=basic_columns))
    cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve/dr9_lrg_{}_1.0.0_photom.fits'.format(field), columns=photom_columns))
    pz = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve/dr9_lrg_{}_1.0.0_pz_new.fits'.format(field)))
    lrgmask = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve/dr9_lrg_{}_1.0.0_lrgmask_v1.fits'.format(field)))
    cat = hstack([cat0, cat, pz, lrgmask], join_type='exact')
    print(len(cat))

    mask = apply_mask(cat, min_nobs, maskbits)
    cat = cat[mask]

    cat['pz_bin'] = np.int16(-1)

    # with warnings.catch_warnings():
    #     warnings.simplefilter("ignore")
    #     gmag = 22.5 - 2.5*np.log10(cat['FLUX_G']/cat['MW_TRANSMISSION_G'])
    #     rmag = 22.5 - 2.5*np.log10(cat['FLUX_R']/cat['MW_TRANSMISSION_R'])
    #     zmag = 22.5 - 2.5*np.log10(cat['FLUX_Z']/cat['MW_TRANSMISSION_Z'])
    #     w1mag = 22.5 - 2.5*np.log10(cat['FLUX_W1']/cat['MW_TRANSMISSION_W1'])
    #     zfibermag = 22.5 - 2.5*np.log10(cat['FIBERFLUX_Z']/cat['MW_TRANSMISSION_Z'])

    ############################## photo-z bins ##############################

    if field=='south':
        pz_cuts = [0.400, 0.540, 0.713, 0.860, 1.020]
    else:
        pz_cuts = [0.400, 0.545, 0.719, 0.851, 1.024]

    for bin_index in range(len(pz_cuts)-1):

        pz_min, pz_max = pz_cuts[bin_index], pz_cuts[bin_index+1]
        mask = (cat['Z_PHOT_MEDIAN']>=pz_min) & (cat['Z_PHOT_MEDIAN']<pz_max)
        cat['pz_bin'][mask] = bin_index+1

    cat_stack.append(cat)

cat = vstack(cat_stack)
# cat = cat[columns_to_keep]
print(len(cat))

output_path = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/main_lrg_minobs_{}_20210913.fits'.format(min_nobs)
print(output_path)
cat.write(output_path, overwrite=False)

