# Example script that load the randoms catalog and applies cuts to match the survey geometry of the mofieid LRG catalog that use DR10 photo-z's
# It matches the quality & footprint cuts as in quality_cuts_dr10_pz.py

import numpy as np
from astropy.table import Table, hstack
import fitsio
import healpy as hp

randoms_index_str = '1-0'

###################### Create the same survey geometry as the original LRG catalogs ######################

# randoms catalog
columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'EBV', 'PHOTSYS']
randoms = Table(fitsio.read(f'/dvs_ro/cfs/cdirs/desi/public/ets/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-{randoms_index_str}.fits', columns=columns))
# LRG mask for the randoms
lrgmask = Table(fitsio.read(f'/dvs_ro/cfs/cdirs/desi/public/papers/c3/lrg_xcorr_2023/v1/catalogs/lrgmask_v1.1/randoms-{randoms_index_str}-lrgmask_v1.1.fits.gz'))
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# DR10 NEXP values row-matched to DR9 randoms
dr10nexp = Table(fitsio.read(f'/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr10_photoz/randoms/randoms-{randoms_index_str}-dr10_nexp.fits.gz'))
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
randoms = hstack([randoms, lrgmask, dr10nexp])
print(len(randoms))

target_min_nobs = 1
target_maskbits = sorted([1, 12, 13])

# DR9 NOBS cut
mask = (randoms['NOBS_G']>=target_min_nobs) & (randoms['NOBS_R']>=target_min_nobs) & (randoms['NOBS_Z']>=target_min_nobs)
randoms = randoms[mask]

mask = np.ones(len(randoms), dtype=bool)
for bit in target_maskbits:
    mask &= (randoms['MASKBITS'] & 2**bit)==0
randoms = randoms[mask]

###################### Create the same survey geometry as the modified LRG catalogs with DR10 photo-z's ######################

min_nobs = 1
max_ebv = 0.15
max_stardens = 2500

# Remove "islands" in the NGC
mask = ~((randoms['DEC']<-10.5) & (randoms['RA']>120) & (randoms['RA']<260))
print('Remove islands', np.sum(mask), np.sum(~mask), np.sum(mask)/len(mask))
randoms = randoms[mask]

# NOBS cut
mask = (randoms['NOBS_G']>=min_nobs) & (randoms['NOBS_R']>=min_nobs) & (randoms['NOBS_Z']>=min_nobs)
print('NOBS', np.sum(mask), np.sum(~mask), np.sum(mask)/len(mask))
randoms = randoms[mask]

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# DR10 NOBS cut
mask = (randoms['DR10_NOBS_G']>=min_nobs) & (randoms['DR10_NOBS_R']>=min_nobs) & (randoms['DR10_NOBS_I']>=min_nobs) & (randoms['DR10_NOBS_Z']>=min_nobs)
print('DR10 NOBS', np.sum(mask), np.sum(mask)/len(mask))
randoms = randoms[mask]

# PHOTSYS
mask = randoms['PHOTSYS']=='S'
print('PHOTSYS=S', np.sum(~mask))
randoms = randoms[mask]

# Remove region with large fraction of bad photo-z's
bad_exposure_ra, bad_exposure_dec = 33.64, 13.38
mask_bad = (randoms['RA']>bad_exposure_ra-1.5) & (randoms['RA']<bad_exposure_ra+1.5)
mask_bad &= (randoms['DEC']>bad_exposure_dec-1.5) & (randoms['DEC']<bad_exposure_dec+1.5)
print('Bad exposure region', np.sum(mask_bad), np.sum(mask_bad)/len(mask_bad))
randoms = randoms[~mask_bad]

# Remove region near LMC where data and randoms do not match
mask_bad = (randoms['RA']>68) & (randoms['RA']<87)
mask_bad &= (randoms['DEC']<-62)
print('Bad LMC region', np.sum(mask_bad), np.sum(mask_bad)/len(mask_bad))
randoms = randoms[~mask_bad]

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Apply LRG mask
mask = randoms['lrg_mask']==0
print('LRG mask', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
randoms = randoms[mask]

# EBV cut
mask = randoms['EBV']<max_ebv
print('EBV', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
randoms = randoms[mask]

# Stellar density cut
stardens = Table(fitsio.read('/dvs_ro/cfs/cdirs/desi/public/papers/c3/lrg_xcorr_2023/v1/misc/pixweight-dr7.1-0.22.0_stardens_64_ring.fits'))  # Stellar density map
stardens_nside = 64
mask = stardens['STARDENS']>=max_stardens
bad_hp_idx = stardens['HPXPIXEL'][mask]
cat_hp_idx = hp.pixelfunc.ang2pix(stardens_nside, randoms['RA'], randoms['DEC'], lonlat=True, nest=False)
mask_bad = np.in1d(cat_hp_idx, bad_hp_idx)
print('STARDENS', np.sum(~mask_bad), np.sum(mask_bad), np.sum(mask_bad)/len(mask_bad))
randoms = randoms[~mask_bad]

print(len(randoms))

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Identify and remove islands
nside = 256
randoms['pix'] = hp.ang2pix(nside, randoms['RA'], randoms['DEC'], nest=False, lonlat=True)
pix_to_keep = np.load('/global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr10_photoz/misc/not_islands_pix_{}.npy'.format(nside))
mask = np.in1d(randoms['pix'], pix_to_keep)
print('Remove islands', np.sum(~mask), np.sum(mask), np.sum(~mask)/len(mask))
randoms = randoms[mask]
print(len(randoms))

