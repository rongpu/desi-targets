# THIS SCRIPT IS OF NO USE SINCE ALL THE NECESSARY COLUMNS ARE ALREADY IN lrg_magnification_{north,south}.fits
# "Parent" catalog with all necessary columns for all objects in the magnified photo-z catalogs

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from desitarget.targets import decode_targetid, encode_targetid


field = 'south'

columns = ['RELEASE', 'BRICKID', 'OBJID', 'RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 
           'MASKBITS', 'EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FIBERFLUX_Z',
           'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1']


targetid_all = []

for magnification in [0.99, 1., 1.01]:
    cat_path = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/lrg_magnification_pz_{}_{:g}.fits'.format(field, magnification)
    cat = Table(fitsio.read(cat_path, columns=['RELEASE', 'BRICKID', 'OBJID']))
    targetid_all.append(encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE']))

targetid_all = np.unique(np.concatenate(targetid_all))

sweep_fns = glob.glob('/global/project/projectdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0/*.fits'.format(field))

sweep_all = []

for sweep_fn in sweep_fns:

    tmp = Table(fitsio.read(sweep_fn, columns=['RELEASE', 'BRICKID', 'OBJID']))
    targetid = encode_targetid(tmp['OBJID'], tmp['BRICKID'], tmp['RELEASE'])
    idx = np.where(np.in1d(targetid, targetid_all))[0]
    targetid = targetid[idx]
    sweep = Table(fitsio.read(sweep_fn, columns=columns, rows=idx))
    sweep['TARGETID'] = targetid

    sweep_all.append(sweep)

sweep_all = vstack(sweep_all)

sweep_all.write('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/lrg_magnification_parent_{}.fits'.format(field))
