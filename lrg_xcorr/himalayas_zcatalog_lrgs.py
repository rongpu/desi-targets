from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

target_bit = 0  # LRG
fn = '/global/cfs/cdirs/desi/spectro/redux/himalayas/zcatalog/ztile-main-dark-cumulative.fits'
cat = Table(fitsio.read(fn, columns=['DESI_TARGET']))
idx = np.where(cat['DESI_TARGET'] & 2**target_bit > 0)[0]
cat = Table(fitsio.read(fn, rows=idx))
print(len(cat), len(np.unique(cat['TARGETID'])))

cat.write('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/himalayas_zcatalog/ztile-main-dark-cumulative-lrg.fits')
