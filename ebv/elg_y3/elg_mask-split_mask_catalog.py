import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

allcat = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/misc/desi_targets_with_desi_ebv_y3/elgmask/all_combined.fits'))
elgmask = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/misc/desi_targets_with_desi_ebv_y3/elgmask/all_combined_elgmask_v1.fits.gz'))
assert len(allcat)==len(elgmask)
allcat = hstack([allcat, elgmask], join_type='exact')

fns = glob.glob('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/misc/desi_targets_with_desi_ebv_y3/*.fits')
fns.sort()

for fn in fns:
    cat = Table(fitsio.read(fn, columns=['TARGETID', 'RA', 'DEC']))
    print(fn, len(cat))
    mask = np.in1d(allcat['TARGETID'], cat['TARGETID'])
    cat1 = allcat[mask].copy()

    # Matching cat1 to cat
    if len(cat)!=len(cat1) or not np.all(np.unique(cat['TARGETID'])==np.unique(cat1['TARGETID'])):
        raise ValueError('cat and cat1 have different TARGETID list')
    t1_reverse_sort = np.array(cat['TARGETID']).argsort().argsort()
    cat1 = cat1[np.argsort(cat1['TARGETID'])[t1_reverse_sort]]
    assert np.all(cat1['TARGETID']==cat['TARGETID'])

    cat1 = cat1[['elg_mask']]
    cat1.write(os.path.join(os.path.dirname(fn), 'elgmask', os.path.basename(fn).replace('.fits', '_elgmask_v1.fits.gz')))

