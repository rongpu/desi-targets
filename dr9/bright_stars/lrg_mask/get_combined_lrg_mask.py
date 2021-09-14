from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits


unwise_maskbits = [0, 1, 2, 3, 4, 6, 7]  # all except the HALO bit
maskbits = [1, 12, 13]  # DESI targeting mask bits

new_mask_dir = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/dev'

randoms_output_dir = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/lrgmask_v1/randoms'

########################################## LRG catalog ##########################################

for field in ['south', 'north']:
    cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve/dr9_lrg_{}_1.0.0_basic.fits'.format(field)))
    cat_wisemask = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve/dr9_lrg_{}_1.0.0_wisemask.fits'.format(field)))
    wmask = np.load(os.path.join(new_mask_dir, 'dr9_lrg_{}_1.0.0_basic-wisemask.npz'.format(field)))
    gmask = np.load(os.path.join(new_mask_dir, 'dr9_lrg_{}_1.0.0_basic-gaiamask.npz'.format(field)))

    cat = hstack([cat, cat_wisemask], join_type='exact')
    cat['wise_mask'] = wmask['wise_mask']
    cat['gaia_mask'] = gmask['gaia_mask']
    cat['gaia_bright_mask'] = gmask['gaia_bright_mask']

    # Apply the targeting maskbits
    mask_ts = np.zeros(len(cat), dtype=bool)
    for bit in maskbits:
        mask_ts |= (cat['MASKBITS'] & 2**bit)>0
    print(np.sum(mask_ts), np.sum(mask_ts)/len(mask_ts))

    # Apply the unWISE maskbits
    mask_unwise = np.zeros(len(cat), dtype=bool)
    for bit in unwise_maskbits:
        mask_unwise |= (cat['WISEMASK_W1'] & 2**bit)>0
    wise_mask = cat['wise_mask'] | mask_unwise
    new_mask = wise_mask | cat['gaia_mask']

    maskbits = [1, 12, 13]

    cat['lrg_mask'] = np.zeros(len(cat), dtype=np.int16)

    # Bit 0: unWISE maskbits
    cat['lrg_mask'][mask_unwise] += 2**0

    # Bit 1: WISE mask
    cat['lrg_mask'][cat['wise_mask']] += 2**1

    # Bit 2: GAIA mask
    cat['lrg_mask'][cat['gaia_mask']] += 2**2

    # Bit 3: GAIA bright mask
    cat['lrg_mask'][cat['gaia_bright_mask']] += 2**3

    # Bit 4: Target selection mask
    cat['lrg_mask'][mask_ts] += 2**4

    cat = cat[['lrg_mask']]
    cat.write('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve/dr9_lrg_{}_1.0.0_lrgmask_v1.fits'.format(field), overwrite=True)

########################################## Randoms ##########################################

randoms_columns = ['MASKBITS', 'PHOTSYS', 'WISEMASK_W1']
randoms_paths = sorted(glob.glob('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-[0-9]*.fits'))

for randoms_path in randoms_paths:

    time_start = time.time()
    
    output_path = os.path.join(randoms_output_dir, os.path.basename(randoms_path).replace('.fits', '-lrgmask_v1.fits'))

    wmask_path = os.path.join(new_mask_dir, os.path.basename(randoms_path).replace('.fits', '-wisemask.npz'))
    gmask_path = os.path.join(new_mask_dir, os.path.basename(randoms_path).replace('.fits', '-gaiamask.npz'))
    
    randoms = Table(fitsio.read(randoms_path, columns=randoms_columns))
    wmask = np.load(wmask_path)
    gmask = np.load(gmask_path)
    
    # Apply the targeting maskbits
    mask_ts = np.zeros(len(randoms), dtype=bool)
    for bit in maskbits:
        mask_ts |= (randoms['MASKBITS'] & 2**bit)>0
    print(np.sum(mask_ts), np.sum(mask_ts)/len(mask_ts))

    # Apply the unWISE maskbits
    mask_unwise = np.zeros(len(randoms), dtype=bool)
    for bit in unwise_maskbits:
        mask_unwise |= (randoms['WISEMASK_W1'] & 2**bit)>0
    wise_mask = wmask['wise_mask'] | mask_unwise
    new_mask = wise_mask | gmask['gaia_mask']

    randoms['lrg_mask'] = np.zeros(len(randoms), dtype=np.int16)

    # Bit 0: unWISE maskbits
    randoms['lrg_mask'][mask_unwise] += 2**0

    # Bit 1: WISE mask
    randoms['lrg_mask'][wmask['wise_mask']] += 2**1

    # Bit 2: GAIA mask
    randoms['lrg_mask'][gmask['gaia_mask']] += 2**2

    # Bit 3: GAIA bright mask
    randoms['lrg_mask'][gmask['gaia_bright_mask']] += 2**3

    # Bit 4: Target selection mask
    randoms['lrg_mask'][mask_ts] += 2**4

    randoms = randoms[['lrg_mask']]    

    randoms.write(output_path, overwrite=True)
    
    print(randoms_path, time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
