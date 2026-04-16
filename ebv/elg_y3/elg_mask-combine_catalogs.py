import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

fns = glob.glob('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/misc/desi_targets_with_desi_ebv_y3/*.fits')
fns.sort()

cat_stack = []
for fn in fns:
    cat = Table(fitsio.read(fn, columns=['TARGETID', 'RA', 'DEC']))
    print(len(cat))
    cat_stack.append(cat)
cat = vstack(cat_stack, join_type='exact')
print(len(cat))


print(len(cat), len(np.unique(cat['TARGETID'])))
_, idx = np.unique(cat['TARGETID'], return_index=True)
cat = cat[idx]
print(len(cat), len(np.unique(cat['TARGETID'])))

cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/misc/desi_targets_with_desi_ebv_y3/elgmask/all_combined.fits')
