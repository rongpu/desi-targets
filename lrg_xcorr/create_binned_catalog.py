from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

columns_to_keep = ['TARGETID', 'RA', 'DEC', 'EBV', 'PIXEL_NOBS_G', 'PIXEL_NOBS_R', 'PIXEL_NOBS_Z', 'MASKBITS', 'PHOTSYS', 'Z_PHOT_MEDIAN', 'lrg_mask', 'pz_bin']

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.1.1/resolve/dr9_lrg_1.1.1_basic.fits'))
photom = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.1.1/resolve/dr9_lrg_1.1.1_photom.fits', columns=['EBV']))
pixel = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.1.1/resolve/dr9_lrg_1.1.1_pixel.fits'))
pz = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.1.1/resolve/dr9_lrg_1.1.1_pz_new.fits', columns=['Z_PHOT_MEDIAN']))
lrgmask = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.1.1/resolve/dr9_lrg_1.1.1_lrgmask_v1.1.fits.gz'))
cat = hstack([cat, photom, pixel, pz, lrgmask], join_type='exact')
print(len(cat))

cat['pz_bin'] = np.int16(-1)

# South
mask0 = cat['PHOTSYS']=='S'
pz_cuts = [0.400, 0.540, 0.713, 0.860, 1.020]
for bin_index in range(len(pz_cuts)-1):
    pz_min, pz_max = pz_cuts[bin_index], pz_cuts[bin_index+1]
    mask = mask0 & (cat['Z_PHOT_MEDIAN']>=pz_min) & (cat['Z_PHOT_MEDIAN']<pz_max)
    cat['pz_bin'][mask] = bin_index+1

# North
mask0 = cat['PHOTSYS']=='N'
pz_cuts = [0.400, 0.545, 0.719, 0.851, 1.024]
for bin_index in range(len(pz_cuts)-1):
    pz_min, pz_max = pz_cuts[bin_index], pz_cuts[bin_index+1]
    mask = mask0 & (cat['Z_PHOT_MEDIAN']>=pz_min) & (cat['Z_PHOT_MEDIAN']<pz_max)
    cat['pz_bin'][mask] = bin_index+1

cat = cat[columns_to_keep]
print(len(cat))

output_path = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/dr9_lrg_1.1.1_pzbins_20221102.fits'
print(output_path)

cat.write(output_path, overwrite=True)
