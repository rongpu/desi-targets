# Subsample v1.0
# More columns for linear weights

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
maskbits = [1, 8, 9, 11, 12, 13]

basic_columns = ['TARGETID', 'RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'PHOTSYS']
photom_columns = ['EBV']
more_columns = ['GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z', 'PSFDEPTH_W1', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z']
columns_to_keep = ['TARGETID', 'EBV', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 'galdepth_gmag_ebv', 'galdepth_rmag_ebv', 'galdepth_zmag_ebv', 'psfdepth_w1mag_ebv']

cat_stack = []
for field in ['north', 'south']:

    cat0 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve/dr9_lrg_{}_1.0.0_basic.fits'.format(field), columns=basic_columns))
    cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve/dr9_lrg_{}_1.0.0_photom.fits'.format(field), columns=photom_columns))
    cat_more = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve/dr9_lrg_{}_1.0.0_more_2.fits'.format(field), columns=more_columns))
    pz = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve/dr9_lrg_{}_1.0.0_pz_new.fits'.format(field)))
    cat = hstack([cat0, cat, pz, cat_more], join_type='exact')
    print(len(cat))

    mask = apply_mask(cat, min_nobs, maskbits)
    cat = cat[mask]

    cat['pz_bin'] = np.int16(-1)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cat['galdepth_gmag_ebv'] = -2.5*(np.log10((5/np.sqrt(cat['GALDEPTH_G'])))-9) - 3.214*cat['EBV']
        cat['galdepth_rmag_ebv'] = -2.5*(np.log10((5/np.sqrt(cat['GALDEPTH_R'])))-9) - 2.165*cat['EBV']
        cat['galdepth_zmag_ebv'] = -2.5*(np.log10((5/np.sqrt(cat['GALDEPTH_Z'])))-9) - 1.211*cat['EBV']
        cat['psfdepth_w1mag_ebv'] = -2.5*(np.log10((5/np.sqrt(cat['PSFDEPTH_W1'])))-9) - 0.184*cat['EBV']
        # cat['psfdepth_w2mag_ebv'] = -2.5*(np.log10((5/np.sqrt(cat['PSFDEPTH_W2'])))-9) - 0.113*cat['EBV']

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
cat = cat[columns_to_keep]
print(len(cat))

output_path = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/main_lrg_minobs_{}_maskbits_{}_20210723_more.fits'.format(min_nobs, ''.join([str(tmp) for tmp in maskbits]))
print(output_path)
cat.write(output_path, overwrite=False)
