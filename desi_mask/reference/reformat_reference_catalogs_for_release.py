from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_lrg_mask_v1.fits'))
cat.rename_columns(['mask_mag', 'radius_south', 'radius_north'], ['MASK_MAG', 'RADIUS_SOUTH', 'RADIUS_NORTH'])
for col in ['MASK_MAG', 'RADIUS_SOUTH', 'RADIUS_NORTH']:
    cat[col] = np.array(cat[col], dtype='float32')
colnames = cat.colnames
data = [np.array(cat[col]) for col in cat.colnames]
fits = fitsio.FITS('/global/cfs/cdirs/desi/users/rongpu/ets/lrg_veto_mask/gaia_lrg_mask_v1.fits.gz', 'rw')
fits.write(data, names=colnames, extname='GAIA_MASK')
fits.close()

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/desi_mask/w1_bright-2mass-lrg_mask_v1.fits'))
cat.rename_columns(['w1ab', 'radius'], ['W1AB', 'RADIUS'])
for col in ['W1MPRO', 'W1AB', 'RADIUS']:
    cat[col] = np.array(cat[col], dtype='float32')
colnames = cat.colnames
data = [np.array(cat[col]) for col in cat.colnames]
fits = fitsio.FITS('/global/cfs/cdirs/desi/users/rongpu/ets/lrg_veto_mask/w1_bright-2mass-lrg_mask_v1.fits.gz', 'rw')
fits.write(data, names=colnames, extname='WISE_MASK')
fits.close()
