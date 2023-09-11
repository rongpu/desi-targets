# Example script that applies the same quality/footprint cuts as in White et al. 2022

import numpy as np
from astropy.table import Table
import fitsio
import healpy as hp

min_nobs = 2
max_ebv = 0.15
max_stardens = 2500

cat = Table(fitsio.read('catalogs/dr9_lrg_pzbins_20230509.fits'))
print(len(cat))

# Remove "islands" in the NGC
mask = ~((cat['DEC']<-10.5) & (cat['RA']>120) & (cat['RA']<260))
print('Remove islands', np.sum(mask), np.sum(~mask), np.sum(mask)/len(mask))
cat = cat[mask]

# NOBS cut
mask = (cat['PIXEL_NOBS_G']>=min_nobs) & (cat['PIXEL_NOBS_R']>=min_nobs) & (cat['PIXEL_NOBS_Z']>=min_nobs)
print('NOBS', np.sum(mask), np.sum(~mask), np.sum(mask)/len(mask))
cat = cat[mask]

# Apply LRG mask
mask = cat['lrg_mask']==0
print('LRG mask', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
cat = cat[mask]

# EBV cut
mask = cat['EBV']<max_ebv
print('EBV', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
cat = cat[mask]

# STARDENS cut
stardens = np.load('misc/pixweight-dr7.1-0.22.0_stardens_64_ring.npy')  # Stellar density map
stardens_nside = 64
mask = stardens>=max_stardens
bad_hp_idx = np.arange(len(stardens))[mask]
cat_hp_idx = hp.pixelfunc.ang2pix(stardens_nside, cat['RA'], cat['DEC'], lonlat=True, nest=False)
mask_bad = np.in1d(cat_hp_idx, bad_hp_idx)
print('STARDENS', np.sum(~mask_bad), np.sum(mask_bad), np.sum(mask_bad)/len(mask_bad))
cat = cat[~mask_bad]

print(len(cat))
