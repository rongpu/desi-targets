from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

# input_dir = '/dvs_ro/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'
input_dir = '/dvs_ro/cfs/cdirs/desi/target/catalogs/dr11/5.1.0/randoms/resolve/'

# n_split = 4
n_split = 1

fns_all = glob.glob(input_dir + 'randoms-[0-9]*')
fns_all.sort()
fns_all = [os.path.basename(tmp) for tmp in fns_all]
print(len(fns_all))

for ii in range(n_split):
    cat = Table(fitsio.read('/pscratch/sd/r/rongpu/tmp/dr11_randoms_masking/randoms_group_{}.fits'.format(ii), columns=['file_index']))
    mm = Table(fitsio.read('/pscratch/sd/r/rongpu/tmp/dr11_randoms_masking/randoms_group_{}-lrgmask_v1.1.fits.gz'.format(ii)))
    if len(cat)!=len(mm):
        raise ValueError
    file_indices = np.unique(cat['file_index'])
    for file_index in file_indices:
        fn = os.path.basename(fns_all[file_index]).replace('.fits', '-lrgmask_v1.1.fits.gz')
        mask = cat['file_index']==file_index
        tmp = mm[mask].copy()
        tmp.write('/global/cfs/cdirs/desi/users/rongpu/desi_mask/randoms_dr11/lrgmask_v1.1/'+fn, overwrite=False)

for ii in range(n_split):
    cat = Table(fitsio.read('/pscratch/sd/r/rongpu/tmp/dr11_randoms_masking/randoms_group_{}.fits'.format(ii), columns=['file_index']))
    mm = Table(fitsio.read('/pscratch/sd/r/rongpu/tmp/dr11_randoms_masking/randoms_group_{}-elgmask_v1.fits.gz'.format(ii)))
    if len(cat)!=len(mm):
        raise ValueError
    file_indices = np.unique(cat['file_index'])
    for file_index in file_indices:
        fn = os.path.basename(fns_all[file_index]).replace('.fits', '-elgmask_v1.fits.gz')
        mask = cat['file_index']==file_index
        tmp = mm[mask].copy()
        tmp.write('/global/cfs/cdirs/desi/users/rongpu/desi_mask/randoms_dr11/elgmask_v1/'+fn, overwrite=False)
