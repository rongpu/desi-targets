# Example script that applies quality cuts on the modified catalogs that use DR10 photo-z's

import numpy as np
from astropy.table import Table
import fitsio
import healpy as hp

min_nobs = 1
max_ebv = 0.15
max_stardens = 2500

cat = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr10_photoz/dr9_lrg_pzbins-dr10_pz.fits'))
# Extended LRG catalog:
# cat = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr10_photoz/dr9_extended_lrg_pzbins-dr10_pz.fits'))
print(len(cat))

# Remove "islands" in the NGC
mask = ~((cat['DEC']<-10.5) & (cat['RA']>120) & (cat['RA']<260))
print('Remove islands', np.sum(mask), np.sum(~mask), np.sum(mask)/len(mask))
cat = cat[mask]

# NOBS cut
mask = (cat['PIXEL_NOBS_G']>=min_nobs) & (cat['PIXEL_NOBS_R']>=min_nobs) & (cat['PIXEL_NOBS_Z']>=min_nobs)
print('NOBS', np.sum(mask), np.sum(~mask), np.sum(mask)/len(mask))
cat = cat[mask]

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# DR10 NOBS cut
mask = (cat['DR10_PIXEL_NOBS_G']>=min_nobs) & (cat['DR10_PIXEL_NOBS_R']>=min_nobs) & (cat['DR10_PIXEL_NOBS_I']>=min_nobs) & (cat['DR10_PIXEL_NOBS_Z']>=min_nobs)
print('DR10 NOBS', np.sum(mask), np.sum(~mask), np.sum(mask)/len(mask))
cat = cat[mask]

# Remove region with large fraction of bad photo-z's
bad_exposure_ra, bad_exposure_dec = 33.64, 13.38
mask_bad = (cat['RA']>bad_exposure_ra-1.5) & (cat['RA']<bad_exposure_ra+1.5)
mask_bad &= (cat['DEC']>bad_exposure_dec-1.5) & (cat['DEC']<bad_exposure_dec+1.5)
print('Bad exposure region', np.sum(mask_bad), np.sum(mask_bad)/len(mask_bad))
cat = cat[~mask_bad]

# Remove region near LMC where data and randoms do not match
mask_bad = (cat['RA']>68) & (cat['RA']<87)
mask_bad &= (cat['DEC']<-62)
print('Bad LMC region', np.sum(mask_bad), np.sum(mask_bad)/len(mask_bad))
cat = cat[~mask_bad]

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Apply LRG mask
mask = cat['lrg_mask']==0
print('LRG mask', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
cat = cat[mask]

# EBV cut
mask = cat['EBV']<max_ebv
print('EBV', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
cat = cat[mask]

# Stellar density cut
stardens = Table(fitsio.read('/dvs_ro/cfs/cdirs/desi/public/papers/c3/lrg_xcorr_2023/v1/misc/pixweight-dr7.1-0.22.0_stardens_64_ring.fits'))  # Stellar density map
stardens_nside = 64
mask = stardens['STARDENS']>=max_stardens
bad_hp_idx = stardens['HPXPIXEL'][mask]
cat_hp_idx = hp.pixelfunc.ang2pix(stardens_nside, cat['RA'], cat['DEC'], lonlat=True, nest=False)
mask_bad = np.in1d(cat_hp_idx, bad_hp_idx)
print('STARDENS', np.sum(~mask_bad), np.sum(mask_bad), np.sum(mask_bad)/len(mask_bad))
cat = cat[~mask_bad]

print(len(cat))

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Identify and remove islands
nside = 256
cat['pix'] = hp.ang2pix(nside, cat['RA'], cat['DEC'], nest=False, lonlat=True)
pix_to_keep = np.load('/global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr10_photoz/misc/not_islands_pix_{}.npy'.format(nside))
mask = np.in1d(cat['pix'], pix_to_keep)
print('Remove islands', np.sum(~mask), np.sum(mask), np.sum(~mask)/len(mask))
cat = cat[mask]
print(len(cat))

