from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

fns = glob.glob('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-[0-9]*')
fns.sort()
fns = [os.path.basename(tmp) for tmp in fns]
print(len(fns))

fns_split, file_indices_split = np.array_split(fns, 4), np.array_split(np.arange(len(fns)), 4)

for ii, fns in enumerate(fns_split):
    cat = Table(fitsio.read('/pscratch/sd/r/rongpu/tmp/randoms_group_{}.fits'.format(ii), columns=['file_index']))
    lrgmask = Table(fitsio.read('/pscratch/sd/r/rongpu/tmp/randoms_group_{}-lrgmask_v1.1.fits.gz'.format(ii)))
    if len(cat)!=len(lrgmask):
        raise ValueError
    file_indices = file_indices_split[ii]
    for file_index in file_indices:
        fn = os.path.basename(fns[file_index]).replace('.fits', '-lrgmask_v1.1.fits.gz')
        mask = cat['file_index']==file_index
        tmp = lrgmask[mask].copy()
        tmp.write('/global/cfs/cdirs/desi/users/rongpu/desi_mask/randoms/lrgmask_v1.1/'+fn)
