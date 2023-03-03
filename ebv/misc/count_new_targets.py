from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp


def apply_mask(cat, min_nobs, maskbits, custom_mask_name):

    # mask = (cat['PIXEL_NOBS_G']>=min_nobs) & (cat['PIXEL_NOBS_R']>=min_nobs) & (cat['PIXEL_NOBS_Z']>=min_nobs)
    mask = np.full(len(cat), True)

    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    # print(np.sum(~mask_clean)/len(mask_clean))

    if custom_mask_name!='':
        mask_col = custom_mask_name[: custom_mask_name.find("mask")]+'_mask'
        mask_clean &= cat[mask_col]==0

    mask &= mask_clean

    return mask


# target_class = 'LRG'
target_class = 'ELG_LOP'

min_nobs = 1
maskbits_dict = {'LRG': [], 'ELG': [], 'ELG_LOP': [], 'QSO': [1, 8, 9, 11, 12, 13], 'BGS_ANY': [1, 13], 'BGS_BRIGHT': [1, 13]}
custom_mask_dict = {'LRG': 'lrgmask_v1.1', 'ELG': 'elgmask_v1', 'ELG_LOP': 'elgmask_v1', 'QSO': '', 'BGS_ANY': '', 'BGS_BRIGHT': ''}
maskbits = maskbits_dict[target_class]
custom_mask_name = custom_mask_dict[target_class]
mask_str = ''.join([str(tmp) for tmp in maskbits])
if custom_mask_name!='':
    mask_str += '_' + custom_mask_name

if target_class=='LRG':
    fn = '/pscratch/sd/r/rongpu/ebv/lrg_targets_hi_ebv.fits'
    maskfn = '/pscratch/sd/r/rongpu/ebv/lrg_targets_hi_ebv_lrgmask_v1.1.fits.gz'
    nexpfn = '/pscratch/sd/r/rongpu/ebv/lrg_targets_hi_ebv_nexp.fits'
elif target_class=='ELG_LOP':
    fn = '/pscratch/sd/r/rongpu/ebv/elg_targets_hi_ebv.fits'
    maskfn = '/pscratch/sd/r/rongpu/ebv/elg_targets_hi_ebv_elgmask_v1.fits.gz'
    nexpfn = '/pscratch/sd/r/rongpu/ebv/elg_targets_hi_ebv_nexp.fits'

cat = Table(fitsio.read(fn))
vetomask = Table(fitsio.read(maskfn))
# nexp = Table(fitsio.read(nexpfn))
# cat = hstack([cat, vetomask, nexp])
cat = hstack([cat, vetomask])

if target_class=='ELG_LOP':
    mask = cat['elglop'].copy()
    print(np.sum(mask)/len(mask))
    cat = cat[mask]

mask = apply_mask(cat, min_nobs, maskbits, custom_mask_name)
print(np.sum(~mask)/len(mask))
cat = cat[mask]

nside = 1024
npix = hp.nside2npix(nside)
pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'], cat['DEC'], lonlat=True)
pix_unique, pix_count = np.unique(pix_allobj, return_counts=True)
hp_table = Table()
hp_table['HPXPIXEL'] = pix_unique
hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_unique, nest=False, lonlat=True)
hp_table['n_obj'] = pix_count
hp_table.write('/pscratch/sd/r/rongpu/ebv/count_map_{}_{}.fits'.format(target_class.lower(), nside))

