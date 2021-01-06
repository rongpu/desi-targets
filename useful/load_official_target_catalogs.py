from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack
import fitsio
# from astropy.io import fits

target_dir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.47.0/targets/main/resolve/dark'
randoms_path = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.47.0/randoms/resolve/randoms-1-0.fits'

target_class = 'LRG'
field = 'south'

if field=='south':
    photsys = 'S'
else:
    photsys = 'N'

target_bits = {'LRG': 0, 'ELG': 1, 'QSO': 2, 'BGS_ANY': 60}
target_bit = target_bits[target_class]

randoms_density = fitsio.read_header(randoms_path, ext=1)['DENSITY']  # randoms per sq. deg.
min_nobs = 1

target_columns = ['BRICKID', 'BRICK_OBJID', 'MORPHTYPE', 'RA', 'DEC', 'EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'FLUX_W1', 'FLUX_W2', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2', 'WISEMASK_W1', 'WISEMASK_W2', 'MASKBITS', 'PHOTSYS', 'DESI_TARGET']
randoms_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'PHOTSYS']


# Load targets
target_path_list = glob.glob(os.path.join(target_dir, 'targets-dark-hp-*.fits'))
cat = []
for target_path in target_path_list:
    print(target_path)
    tmp = fitsio.read(target_path, columns=['DESI_TARGET', 'PHOTSYS'])
    mask = ((tmp["DESI_TARGET"] & (2**target_bit))!=0) & (tmp['PHOTSYS']==photsys)
    idx = np.where(mask)[0]
    if len(idx)==0:
        continue
    print(len(idx)/len(tmp), len(idx), len(tmp))
    cat.append(Table(fitsio.read(target_path, columns=target_columns, rows=idx)))
cat = vstack(cat)

# Apply WISE mask
# Mask bits already applied: [1, 12, 13]
maskbits = [1, 8, 9, 12, 13]
mask_clean = np.ones(len(cat), dtype=bool)
for bit in maskbits:
    mask_clean &= (cat['MASKBITS'] & 2**bit)==0
print(np.sum(~mask_clean), np.sum(~mask_clean)/len(mask_clean))
cat = cat[mask_clean]
print(len(cat))

mask = (cat['NOBS_G']>=min_nobs) & (cat['NOBS_R']>=min_nobs) & (cat['NOBS_Z']>=min_nobs)
cat = cat[mask]
print(len(cat))
print()


# Load randoms
randoms = Table(fitsio.read(randoms_path, columns=randoms_columns))
print(len(randoms))

mask = (randoms['PHOTSYS']==photsys)
randoms = randoms[mask]

mask = (randoms['NOBS_G']>=min_nobs) & (randoms['NOBS_R']>=min_nobs) & (randoms['NOBS_Z']>=min_nobs)
randoms = randoms[mask]
print(len(randoms))

# Apply masks
maskbits = [1, 8, 9, 12, 13]
mask_clean = np.ones(len(randoms), dtype=bool)
for bit in maskbits:
    mask_clean &= (randoms['MASKBITS'] & 2**bit)==0
print(np.sum(~mask_clean), np.sum(~mask_clean)/len(mask_clean))
randoms = randoms[mask_clean]

print('Area = {:.1f} sq. deg.'.format(len(randoms)/randoms_density))
print('Target density = {:.1f} per sq. deg.'.format(len(cat)/(len(randoms)/randoms_density)))
