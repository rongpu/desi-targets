# Create the catalog for Boryana's kSZ project
# First need to get the DR10 pixel NOBS for the LRG catalog
# python read_pixel_nexp.py --dr 10 --input /global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr9_lrg_pzbins_20230509.fits --output /global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr10_photoz/dr9_lrg_dr10_pixel.fits.gz
# python read_pixel_nexp.py --dr 10 --input /global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr9_extended_lrg_pzbins_20230509.fits --output /global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr10_photoz/dr9_extended_lrg_dr10_pixel.fits.gz

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

############################################# Main LRGs #############################################

dr9 = Table(fitsio.read('/global/cfs/cdirs/desi/target/analysis/truth/dr10.1/south/matched/dr9_lrg_pzbins-match.fits'))
dr10 = Table(fitsio.read('/global/cfs/cdirs/desi/target/analysis/truth/dr10.1/south/matched/ls-dr10.1-dr9_lrg_pzbins-match.fits',
                         columns=['Z_PHOT_MEAN_I', 'Z_PHOT_MEDIAN_I', 'Z_PHOT_STD_I', 'Z_PHOT_L68_I', 'Z_PHOT_U68_I', 'Z_PHOT_L95_I', 'Z_PHOT_U95_I']))
dr10.rename_columns(['Z_PHOT_MEAN_I', 'Z_PHOT_MEDIAN_I', 'Z_PHOT_STD_I', 'Z_PHOT_L68_I', 'Z_PHOT_U68_I', 'Z_PHOT_L95_I', 'Z_PHOT_U95_I'],
                    ['DR10_Z_PHOT_MEAN_I', 'DR10_Z_PHOT_MEDIAN_I', 'DR10_Z_PHOT_STD_I', 'DR10_Z_PHOT_L68_I', 'DR10_Z_PHOT_U68_I', 'DR10_Z_PHOT_L95_I', 'DR10_Z_PHOT_U95_I'])

dr10_pixel = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr10_photoz/dr9_lrg_dr10_pixel.fits.gz'))
matched = np.load('/global/cfs/cdirs/desi/target/analysis/truth/dr10.1/south/allobjects/ls-dr10.1-dr9_lrg_pzbins.npy')
dr10_pixel = dr10_pixel[matched]
dr10_pixel.rename_columns(['PIXEL_NOBS_G', 'PIXEL_NOBS_R', 'PIXEL_NOBS_I', 'PIXEL_NOBS_Z'], ['DR10_PIXEL_NOBS_G', 'DR10_PIXEL_NOBS_R', 'DR10_PIXEL_NOBS_I', 'DR10_PIXEL_NOBS_Z'])

cat = dr9.copy()
cat = hstack([cat, dr10_pixel])
cat = hstack([cat, dr10])
print(len(cat))

# Exclude objects with Northern photometry
mask = cat['PHOTSYS']=='S'
print(np.sum(mask)/len(mask))
cat = cat[mask]
print(len(cat))

cat['DR10_pz_bin'] = np.int16(-1)

pz_cuts_south = [0.400, 0.540, 0.713, 0.860, 1.020]
for bin_index in range(len(pz_cuts_south)-1):
    pz_min, pz_max = pz_cuts_south[bin_index], pz_cuts_south[bin_index+1]
    mask = (cat['DR10_Z_PHOT_MEDIAN_I']>=pz_min) & (cat['DR10_Z_PHOT_MEDIAN_I']<pz_max)
    cat['DR10_pz_bin'][mask] = bin_index+1

cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr10_photoz/dr9_lrg_pzbins-dr10_pz.fits', overwrite=True)

############################################# Extended LRGs #############################################

dr9 = Table(fitsio.read('/global/cfs/cdirs/desi/target/analysis/truth/dr10.1/south/matched/dr9_extended_lrg_pzbins-match.fits'))
dr10 = Table(fitsio.read('/global/cfs/cdirs/desi/target/analysis/truth/dr10.1/south/matched/ls-dr10.1-dr9_extended_lrg_pzbins-match.fits',
                         columns=['Z_PHOT_MEAN_I', 'Z_PHOT_MEDIAN_I', 'Z_PHOT_STD_I', 'Z_PHOT_L68_I', 'Z_PHOT_U68_I', 'Z_PHOT_L95_I', 'Z_PHOT_U95_I']))
dr10.rename_columns(['Z_PHOT_MEAN_I', 'Z_PHOT_MEDIAN_I', 'Z_PHOT_STD_I', 'Z_PHOT_L68_I', 'Z_PHOT_U68_I', 'Z_PHOT_L95_I', 'Z_PHOT_U95_I'],
                    ['DR10_Z_PHOT_MEAN_I', 'DR10_Z_PHOT_MEDIAN_I', 'DR10_Z_PHOT_STD_I', 'DR10_Z_PHOT_L68_I', 'DR10_Z_PHOT_U68_I', 'DR10_Z_PHOT_L95_I', 'DR10_Z_PHOT_U95_I'])

dr10_pixel = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr10_photoz/dr9_extended_lrg_dr10_pixel.fits.gz'))
matched = np.load('/global/cfs/cdirs/desi/target/analysis/truth/dr10.1/south/allobjects/ls-dr10.1-dr9_extended_lrg_pzbins.npy')
dr10_pixel = dr10_pixel[matched]
dr10_pixel.rename_columns(['PIXEL_NOBS_G', 'PIXEL_NOBS_R', 'PIXEL_NOBS_I', 'PIXEL_NOBS_Z'], ['DR10_PIXEL_NOBS_G', 'DR10_PIXEL_NOBS_R', 'DR10_PIXEL_NOBS_I', 'DR10_PIXEL_NOBS_Z'])

cat = dr9.copy()
cat = hstack([cat, dr10_pixel])
cat = hstack([cat, dr10])
print(len(cat))

# Exclude objects with Northern photometry
mask = cat['PHOTSYS']=='S'
print(np.sum(mask)/len(mask))
cat = cat[mask]
print(len(cat))

cat['DR10_pz_bin'] = np.int16(-1)

pz_cuts_south = [0.400, 0.540, 0.713, 0.860, 1.00]
for bin_index in range(len(pz_cuts_south)-1):
    pz_min, pz_max = pz_cuts_south[bin_index], pz_cuts_south[bin_index+1]
    mask = (cat['DR10_Z_PHOT_MEDIAN_I']>=pz_min) & (cat['DR10_Z_PHOT_MEDIAN_I']<pz_max)
    cat['DR10_pz_bin'][mask] = bin_index+1

cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr10_photoz/dr9_extended_lrg_pzbins-dr10_pz.fits', overwrite=True)
