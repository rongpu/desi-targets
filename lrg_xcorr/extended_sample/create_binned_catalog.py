from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

columns_to_keep = ['TARGETID', 'RA', 'DEC', 'EBV', 'PIXEL_NOBS_G', 'PIXEL_NOBS_R', 'PIXEL_NOBS_Z', 'MASKBITS', 'PHOTSYS', 'Z_PHOT_MEDIAN', 'lrg_mask', 'pz_bin']

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/more/dr9_extended_lrg_0.49.0_basic.fits'))
pz = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/more/dr9_extended_lrg_0.49.0_pz_new.fits', columns=['Z_PHOT_MEDIAN']))
lrgmask = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/more/dr9_extended_lrg_0.49.0_lrgmask_v1.1.fits.gz'))
photom = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/more/dr9_extended_lrg_0.49.0_photom.fits', columns=['FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_Z', 'EBV']))
pixel = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/more/dr9_extended_lrg_0.49.0_pixel.fits'))
cat = hstack([cat, pz, lrgmask, photom, pixel], join_type='exact')
print(len(cat))

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    cat['gmag'] = 22.5 - 2.5*np.log10(cat['FLUX_G']) - 3.214 * cat['EBV']
    cat['rmag'] = 22.5 - 2.5*np.log10(cat['FLUX_R']) - 2.165 * cat['EBV']
    cat['zmag'] = 22.5 - 2.5*np.log10(cat['FLUX_Z']) - 1.211 * cat['EBV']
    cat['w1mag'] = 22.5 - 2.5*np.log10(cat['FLUX_W1']) - 0.184 * cat['EBV']
    cat['w2mag'] = 22.5 - 2.5*np.log10(cat['FLUX_W2']) - 0.113 * cat['EBV']
    cat['zfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_Z']) - 1.211 * cat['EBV']

# Adjust for North-South photometry offsets
mask_north = cat['PHOTSYS']=='N'
mask_south = cat['PHOTSYS']=='S'
north_sliding_cut = (cat['rmag'] - cat['w1mag'] > (cat['w1mag'] - 17.44) * 1.8)
north_sliding_cut |= ((cat['rmag']-cat['w1mag'])>3.2)
mask_valid = mask_north & north_sliding_cut
mask_valid |= mask_south & (cat['zfibermag']<21.96)
print(np.sum(mask_valid)/len(mask_valid))
print(np.sum(mask_valid & mask_north)/np.sum(mask_north))
print(np.sum(mask_valid & mask_south)/np.sum(mask_south))

cat['pz_bin'] = np.int16(-1)

pz_cuts_south = [0.400, 0.540, 0.713, 0.860, 1.00]
pz_cuts_north = [0.400, 0.545, 0.719, 0.854, 1.01]

# South
mask0 = mask_valid & (cat['PHOTSYS']=='S')
for bin_index in range(len(pz_cuts_south)-1):
    pz_min, pz_max = pz_cuts_south[bin_index], pz_cuts_south[bin_index+1]
    mask = mask0 & (cat['Z_PHOT_MEDIAN']>=pz_min) & (cat['Z_PHOT_MEDIAN']<pz_max)
    cat['pz_bin'][mask] = bin_index+1

# North
mask0 = mask_valid & (cat['PHOTSYS']=='N')
for bin_index in range(len(pz_cuts_north)-1):
    pz_min, pz_max = pz_cuts_north[bin_index], pz_cuts_north[bin_index+1]
    mask = mask0 & (cat['Z_PHOT_MEDIAN']>=pz_min) & (cat['Z_PHOT_MEDIAN']<pz_max)
    cat['pz_bin'][mask] = bin_index+1

cat = cat[columns_to_keep]
print(len(cat))

cat['extended_lrg'] = mask_valid.copy()

output_path = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/dr9_extended_lrg_0.49.0_pzbins_20230120.fits'
print(output_path)

cat.write(output_path, overwrite=True)
