# Compute number count slope of the SV3 LRGs for Mehdie

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio

import yaml

ff_factor = True  # Include for the fiberflux factor

min_nobs = 2
# maskbits = [1, 8, 9, 11, 12, 13]

columns = ['TARGETID', 'TYPE', 'RA', 'DEC', 'EBV',
'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FIBERFLUX_Z',
'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1',
'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1',
'GAIA_PHOT_G_MEAN_MAG', 'FIBERTOTFLUX_Z',
'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS']

counts = dict()

for field in ['north', 'south']:

    cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/extended_lrg_magnification_{}.fits'.format(field), columns=columns))
    cat1 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/extended_lrg_magnification_{}_fiberflux.fits'.format(field)))
    cat2 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/extended_lrg_magnification_{}_pixel.fits'.format(field)))
    cat3 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/extended_lrg_magnification_{}_lrgmask_v1.1.fits.gz'.format(field)))
    cat = hstack([cat, cat1, cat2, cat3])

    if field=='north':
        cat['PHOTSYS'] = 'N'
    else:
        cat['PHOTSYS'] = 'S'

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cat['gmag'] = 22.5 - 2.5*np.log10(cat['FLUX_G']) - 3.214 * cat['EBV']
        cat['rmag'] = 22.5 - 2.5*np.log10(cat['FLUX_R']) - 2.165 * cat['EBV']
        cat['zmag'] = 22.5 - 2.5*np.log10(cat['FLUX_Z']) - 1.211 * cat['EBV']
        cat['w1mag'] = 22.5 - 2.5*np.log10(cat['FLUX_W1']) - 0.184 * cat['EBV']
        cat['zfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_Z']) - 1.211 * cat['EBV']

    if field=='north':
        mask_ns = (cat['DEC']>32.375)
    else:
        mask_ns = ((cat['DEC']<=32.375) | (cat['RA']<104) | (cat['RA']>280)) & (cat['DEC']>-29)
    cat = cat[mask_ns]

    mask_quality = np.full(len(cat), True)

    mask_quality &= (cat['FLUX_IVAR_R'] > 0) & (cat['FLUX_R'] > 0)   # ADM quality in r.
    mask_quality &= (cat['FLUX_IVAR_Z'] > 0) & (cat['FLUX_Z'] > 0) & (cat['FIBERFLUX_Z'] > 0)   # ADM quality in z.
    mask_quality &= (cat['FLUX_IVAR_W1'] > 0) & (cat['FLUX_W1'] > 0)  # ADM quality in W1.

    mask_quality &= (cat['GAIA_PHOT_G_MEAN_MAG'] == 0) | (cat['GAIA_PHOT_G_MEAN_MAG'] > 18)  # remove bright GAIA sources

    # ADM remove stars with zfibertot < 17.5 that are missing from GAIA.
    mask_quality &= cat['FIBERTOTFLUX_Z'] < 10**(-0.4*(17.5-22.5))

    # mask_quality &= (cat['NOBS_G']>=min_nobs) & (cat['NOBS_R']>=min_nobs) & (cat['NOBS_Z']>=min_nobs)
    mask_quality &= (cat['PIXEL_NOBS_G']>=min_nobs) & (cat['PIXEL_NOBS_R']>=min_nobs) & (cat['PIXEL_NOBS_Z']>=min_nobs)

    # Apply masks
    # mask_clean = np.ones(len(cat), dtype=bool)
    # for bit in maskbits:
    #     mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    # print(np.sum(~mask_clean), np.sum(~mask_clean)/len(mask_clean))
    mask_clean = cat['lrg_mask']==0
    print(np.sum(mask_clean)/len(mask_clean))
    mask_quality &= mask_clean

    for magnification in [0.99, 1., 1.01]:

        gmag = 22.5 - 2.5 * np.log10((cat['FLUX_G'] * magnification / cat['MW_TRANSMISSION_G']).clip(1e-7))
        rmag = 22.5 - 2.5 * np.log10((cat['FLUX_R'] * magnification / cat['MW_TRANSMISSION_R']).clip(1e-7))
        zmag = 22.5 - 2.5 * np.log10((cat['FLUX_Z'] * magnification / cat['MW_TRANSMISSION_Z']).clip(1e-7))
        w1mag = 22.5 - 2.5 * np.log10((cat['FLUX_W1'] * magnification / cat['MW_TRANSMISSION_W1']).clip(1e-7))

        if ff_factor:
            zfibermag = 22.5 - 2.5 * np.log10((cat['FIBERFLUX_Z'] * (1 + (magnification-1) * cat['ff_factor']) / cat['MW_TRANSMISSION_Z']).clip(1e-7))
        else:
            zfibermag = 22.5 - 2.5 * np.log10((cat['FIBERFLUX_Z'] * (1 + (magnification-1) * 1               ) / cat['MW_TRANSMISSION_Z']).clip(1e-7))

        mask_lrg = np.full(len(cat), True)

        if field=='south':
            mask_lrg &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
            mask_lrg &= (zfibermag < 21.7)                   # faint limit

            lrg_mask_sliding = rmag - w1mag > (w1mag - 17.26) * 1.8  # sliding IR cut
            lrg_mask_sliding &= rmag - w1mag > (w1mag - 16.36) * 1.  # low-z sliding IR cut
            lrg_mask_sliding |= rmag - w1mag > 3.29
            mask_lrg &= lrg_mask_sliding

            mask_lowz = (gmag - rmag > 1.3) & ( (gmag - rmag) > -1.55 * (rmag - w1mag) + 3.13)
            mask_lowz |= (rmag - w1mag > 1.8)
            mask_lrg &= mask_lowz
        else:
            mask_lrg &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
            mask_lrg &= (zfibermag < 21.72)                   # faint limit

            lrg_mask_sliding = rmag - w1mag > (w1mag - 17.24) * 1.83  # sliding IR cut
            lrg_mask_sliding &= rmag - w1mag > (w1mag - 16.33) * 1.  # low-z sliding IR cut
            lrg_mask_sliding |= rmag - w1mag > 3.39
            mask_lrg &= lrg_mask_sliding

            mask_lowz = (gmag - rmag > 1.34) & ((gmag - rmag) > -1.55 * (rmag - w1mag) + 3.23)
            mask_lowz |= (rmag - w1mag > 1.8)
            mask_lrg &= mask_lowz

        counts['{}_all_{:.3f}'.format(field, magnification)] = int(np.sum(mask_quality & mask_lrg))

fn = 'sv3_lrg_counts'
if not ff_factor:
    fn += '_no_ff'

with open('/global/cfs/cdirs/desi/users/rongpu/lrg_xcorr/magnification/counts/{}.txt'.format(fn), "w") as f:
    yaml.dump(counts, f)
