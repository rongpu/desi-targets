from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

min_nobs = 1

###################################### Randoms ######################################

randoms_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'PHOTSYS']

for index in range(4):
    randoms_path = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-{}.fits'.format(index)
    print(randoms_path)
    randoms = Table(fitsio.read(randoms_path, columns=randoms_columns))
    for custom_mask_name in ['lrgmask_v1.1', 'elgmask_v1']:
        mask_dir = os.path.join('/global/cfs/cdirs/desi/users/rongpu/desi_mask/randoms/', custom_mask_name)
        mask_path = os.path.join(mask_dir, os.path.basename(randoms_path).replace('.fits', '-{}.fits.gz'.format(custom_mask_name)))
        custom_mask = Table(fitsio.read(mask_path))
        randoms = hstack([randoms, custom_mask], join_type='exact')

    mask = (randoms['NOBS_G']>=min_nobs) & (randoms['NOBS_R']>=min_nobs) & (randoms['NOBS_Z']>=min_nobs)
    randoms = randoms[mask]

    # BGS mask
    maskbits = [1, 13]
    mask_clean = np.ones(len(randoms), dtype=bool)
    for bit in maskbits:
        mask_clean &= (randoms['MASKBITS'] & 2**bit)==0
    randoms['bgs_mask'] = mask_clean.copy()

    # LRG mask
    randoms['lrg_mask'] = randoms['lrg_mask']==0

    # ELG mask
    randoms['elg_mask'] = randoms['elg_mask']==0

    # QSO mask
    maskbits = [1, 8, 9, 11, 12, 13]
    mask_clean = np.ones(len(randoms), dtype=bool)
    for bit in maskbits:
        mask_clean &= (randoms['MASKBITS'] & 2**bit)==0
    randoms['qso_mask'] = mask_clean.copy()

    randoms = randoms[['RA', 'DEC', 'bgs_mask', 'lrg_mask', 'elg_mask', 'qso_mask', 'PHOTSYS']]
    randoms.write('/global/cfs/cdirs/cosmo/www/temp/rongpu/misc/desi_targets_and_randoms/randoms-1-{}-desi.fits'.format(index), overwrite=True)

###################################### Targets ######################################

main_dir = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.1.1/resolve'
columns = ['RA', 'DEC', 'PHOTSYS', 'EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z']

########################## BGS ##########################

print('BGS')

maskbits = [1, 13]
tmp = Table(fitsio.read(os.path.join(main_dir, 'dr9_bgs_any_1.1.1_basic.fits')))
tmp1 = Table(fitsio.read(os.path.join(main_dir, 'dr9_bgs_any_1.1.1_photom.fits')))
tmp3 = Table(fitsio.read(os.path.join(main_dir, 'dr9_bgs_any_1.1.1_pixel.fits')))
cat = hstack([tmp, tmp1, tmp3], join_type='exact')
print(len(cat))

cat['BGS_FAINT'] = cat['BGS_TARGET'] & 2**0 > 0
cat['BGS_BRIGHT'] = cat['BGS_TARGET'] & 2**1 > 0

mask = (cat['PIXEL_NOBS_G']>=min_nobs) & (cat['PIXEL_NOBS_R']>=min_nobs) & (cat['PIXEL_NOBS_Z']>=min_nobs)
cat = cat[mask]
print(len(cat))

mask_clean = np.ones(len(cat), dtype=bool)
for bit in maskbits:
    mask_clean &= (cat['MASKBITS'] & 2**bit)==0
cat = cat[mask_clean]
print(len(cat))

cat = cat[columns+['BGS_FAINT', 'BGS_BRIGHT']]

cat.write('/global/cfs/cdirs/cosmo/www/temp/rongpu/misc/desi_targets_and_randoms/bgs.fits')

########################## LRG ##########################

print('LRG')

maskbits = []
custom_mask_name = 'lrgmask_v1.1'
tmp = Table(fitsio.read(os.path.join(main_dir, 'dr9_lrg_1.1.1_basic.fits')))
tmp1 = Table(fitsio.read(os.path.join(main_dir, 'dr9_lrg_1.1.1_photom.fits')))
tmp2 = Table(fitsio.read(os.path.join(main_dir, 'dr9_lrg_1.1.1_lrgmask_v1.1.fits.gz')))
tmp3 = Table(fitsio.read(os.path.join(main_dir, 'dr9_lrg_1.1.1_pixel.fits')))
cat = hstack([tmp, tmp1, tmp2, tmp3], join_type='exact')
print(len(cat))

mask = (cat['PIXEL_NOBS_G']>=min_nobs) & (cat['PIXEL_NOBS_R']>=min_nobs) & (cat['PIXEL_NOBS_Z']>=min_nobs)
cat = cat[mask]
print(len(cat))

mask_clean = np.ones(len(cat), dtype=bool)
for bit in maskbits:
    mask_clean &= (cat['MASKBITS'] & 2**bit)==0
mask_clean &= cat['lrg_mask']==0
cat = cat[mask_clean]
print(len(cat))

cat = cat[columns]

cat.write('/global/cfs/cdirs/cosmo/www/temp/rongpu/misc/desi_targets_and_randoms/lrg.fits')

########################## ELG ##########################

print('ELG')

maskbits = []
custom_mask_name = 'elgmask_v1'
tmp = Table(fitsio.read(os.path.join(main_dir, 'dr9_elg_1.1.1_basic.fits')))
tmp1 = Table(fitsio.read(os.path.join(main_dir, 'dr9_elg_1.1.1_photom.fits')))
tmp2 = Table(fitsio.read(os.path.join(main_dir, 'dr9_elg_1.1.1_elgmask_v1.fits.gz')))
tmp3 = Table(fitsio.read(os.path.join(main_dir, 'dr9_elg_1.1.1_pixel.fits')))
cat = hstack([tmp, tmp1, tmp2, tmp3], join_type='exact')
print(len(cat))

mask = (cat['PIXEL_NOBS_G']>=min_nobs) & (cat['PIXEL_NOBS_R']>=min_nobs) & (cat['PIXEL_NOBS_Z']>=min_nobs)
cat = cat[mask]
print(len(cat))

mask_clean = np.ones(len(cat), dtype=bool)
for bit in maskbits:
    mask_clean &= (cat['MASKBITS'] & 2**bit)==0
mask_clean &= cat['elg_mask']==0
cat = cat[mask_clean]
print(len(cat))

cat = cat[columns]

cat.write('/global/cfs/cdirs/cosmo/www/temp/rongpu/misc/desi_targets_and_randoms/elg.fits')

########################## QSO ##########################

print('QSO')

maskbits = [1, 8, 9, 11, 12, 13]
tmp = Table(fitsio.read(os.path.join(main_dir, 'dr9_qso_1.1.1_basic.fits')))
tmp1 = Table(fitsio.read(os.path.join(main_dir, 'dr9_qso_1.1.1_photom.fits')))
tmp3 = Table(fitsio.read(os.path.join(main_dir, 'dr9_qso_1.1.1_pixel.fits')))
cat = hstack([tmp, tmp1, tmp3], join_type='exact')
print(len(cat))

mask = (cat['PIXEL_NOBS_G']>=min_nobs) & (cat['PIXEL_NOBS_R']>=min_nobs) & (cat['PIXEL_NOBS_Z']>=min_nobs)
cat = cat[mask]
print(len(cat))

mask_clean = np.ones(len(cat), dtype=bool)
for bit in maskbits:
    mask_clean &= (cat['MASKBITS'] & 2**bit)==0
cat = cat[mask_clean]
print(len(cat))

cat = cat[columns]

cat.write('/global/cfs/cdirs/cosmo/www/temp/rongpu/misc/desi_targets_and_randoms/qso.fits')

