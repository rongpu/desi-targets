from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from desispec.io.util import write_bintable

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_lrg_mask_v1.fits'))
cat.rename_columns(['mask_mag', 'radius_south', 'radius_north'], ['MASK_MAG', 'RADIUS_SOUTH', 'RADIUS_NORTH'])
for col in ['MASK_MAG', 'RADIUS_SOUTH', 'RADIUS_NORTH']:
    cat[col] = np.array(cat[col], dtype='float32')
units = {'RA': 'degrees', 'DEC': 'degrees', 'MASK_MAG': 'AB magnitudes', 'RADIUS_SOUTH': 'arcsec', 'RADIUS_NORTH': 'arcsec'}
write_bintable('/global/cfs/cdirs/desi/users/rongpu/ets/lrg_veto_mask/gaia_lrg_mask_v1.fits.gz', cat, extname='GAIA_MASK', units=units, clobber=True)

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/desi_mask/w1_bright-2mass-lrg_mask_v1.fits'))
cat.rename_columns(['w1ab', 'radius'], ['W1AB', 'RADIUS'])
for col in ['W1MPRO', 'W1AB', 'RADIUS']:
    cat[col] = np.array(cat[col], dtype='float32')
units = {'RA': 'degrees', 'DEC': 'degrees', 'W1MPRO': 'Vega magnitudes', 'W1AB': 'AB magnitudes', 'RADIUS': 'arcsec'}
write_bintable('/global/cfs/cdirs/desi/users/rongpu/ets/lrg_veto_mask/w1_bright-2mass-lrg_mask_v1.fits.gz', cat, extname='WISE_MASK', units=units, clobber=True)
