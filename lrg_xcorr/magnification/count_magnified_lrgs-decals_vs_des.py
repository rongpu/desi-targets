from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio

import yaml

pz_magnification = True  # Include the photo-z shift
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

field = 'south'

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/main_lrg_magnification_{}.fits'.format(field), columns=columns))
cat1 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/main_lrg_magnification_{}_fiberflux.fits'.format(field)))
cat2 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/main_lrg_magnification_{}_pixel_nobs.fits'.format(field)))
cat3 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/main_lrg_magnification_{}_lrgmask_v1.1.fits.gz'.format(field)))
cat = hstack([cat, cat1, cat2, cat3])

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
mask_quality &= (cat['NGOOD_G']>=min_nobs) & (cat['NGOOD_R']>=min_nobs) & (cat['NGOOD_Z']>=min_nobs)

# Apply masks
# mask_clean = np.ones(len(cat), dtype=bool)
# for bit in maskbits:
#     mask_clean &= (cat['MASKBITS'] & 2**bit)==0
# print(np.sum(~mask_clean), np.sum(~mask_clean)/len(mask_clean))
mask_clean = cat['lrg_mask']==0
print(np.sum(mask_clean)/len(mask_clean))
mask_quality &= mask_clean

cat = cat[mask_quality]

sys.path.append(os.path.expanduser('~/git/desi-targets/useful'))
from isdes import get_isdes
cat['isdes'] = get_isdes(cat['RA'], cat['DEC'])
cat_all = cat.copy()

for region in ['DECaLS', 'DES']:

    if region=='DECaLS':
        mask_region = ~cat_all['isdes']
    else:
        mask_region = cat_all['isdes'].copy()
    cat = cat_all[mask_region].copy()

    for magnification in [0.99, 1., 1.01]:

        gmag = 22.5 - 2.5 * np.log10((cat['FLUX_G'] * magnification / cat['MW_TRANSMISSION_G']).clip(1e-7))
        rmag = 22.5 - 2.5 * np.log10((cat['FLUX_R'] * magnification / cat['MW_TRANSMISSION_R']).clip(1e-7))
        zmag = 22.5 - 2.5 * np.log10((cat['FLUX_Z'] * magnification / cat['MW_TRANSMISSION_Z']).clip(1e-7))
        w1mag = 22.5 - 2.5 * np.log10((cat['FLUX_W1'] * magnification / cat['MW_TRANSMISSION_W1']).clip(1e-7))

        if ff_factor:
            zfibermag = 22.5 - 2.5 * np.log10((cat['FIBERFLUX_Z'] * (1 + (magnification-1) * cat['ff_factor']) / cat['MW_TRANSMISSION_Z']).clip(1e-7))
        else:
            zfibermag = 22.5 - 2.5 * np.log10((cat['FIBERFLUX_Z'] * (1 + (magnification-1) * 1               ) / cat['MW_TRANSMISSION_Z']).clip(1e-7))

        if pz_magnification:
            pz = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/main_lrg_magnification_pz_{}_{:g}.fits'.format(field, magnification)))
        else:
            pz = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/main_lrg_magnification_pz_{}_1.fits'.format(field)))
        pz = pz[mask_ns]
        pz = pz[mask_quality]
        pz = pz[mask_region]

        if len(pz)!=len(cat):
            raise ValueError

        mask_lrg = np.full(len(cat), True)

        if field=='south':
            mask_lrg &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
            mask_lrg &= zfibermag < 21.6                   # faint limit
            mask_lrg &= (gmag - w1mag > 2.9) | (rmag - w1mag > 1.8)  # low-z cuts
            mask_lrg &= (
                ((rmag - w1mag > (w1mag - 17.14) * 1.8)
                 & (rmag - w1mag > (w1mag - 16.33) * 1.))
                | (rmag - w1mag > 3.3)
            )  # double sliding cuts and high-z extension
        else:
            mask_lrg &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
            mask_lrg &= zfibermag < 21.61                   # faint limit
            mask_lrg &= (gmag - w1mag > 2.97) | (rmag - w1mag > 1.8)  # low-z cuts
            mask_lrg &= (
                ((rmag - w1mag > (w1mag - 17.13) * 1.83)
                 & (rmag - w1mag > (w1mag - 16.31) * 1.))
                | (rmag - w1mag > 3.4)
            )  # double sliding cuts and high-z extension

        counts['{}_all_{:.3f}'.format(region, magnification)] = int(np.sum(mask_lrg))

        if field=='south':
            pz_cuts = [0.400, 0.540, 0.713, 0.860, 1.020]
        else:
            pz_cuts = [0.400, 0.545, 0.719, 0.851, 1.024]

        for bin_index in range(len(pz_cuts)-1):

            pz_min, pz_max = pz_cuts[bin_index], pz_cuts[bin_index+1]
            mask_pz = (pz['Z_PHOT_MEDIAN']>=pz_min) & (pz['Z_PHOT_MEDIAN']<pz_max)

            counts['{}_bin_{}_{:.3f}'.format(region, bin_index+1, magnification)] = int(np.sum(mask_lrg & mask_pz))

fn = 'main_lrg_counts_decals_vs_des'
if not pz_magnification:
    fn += '_no_pz_mag'
if not ff_factor:
    fn += '_no_ff'

with open('/global/cfs/cdirs/desi/users/rongpu/lrg_xcorr/magnification/counts/{}.txt'.format(fn), "w") as f:
    yaml.dump(counts, f)
