from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

# parent = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/misc/cosmos_lowz_20220404/cosmos_lowz_sample_v0.6.fits'))
parent = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/desi2/cosmos_lowz/targets/cosmos_lowz_targets.fits'))
print(len(parent))
mask = parent['is_target'].copy()
parent = parent[mask]
print(len(parent))

parent.rename_columns(['TYPE'], ['MORPHTYPE'])

columns_1 = ['TARGETID', 'CHI2', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'DELTACHI2']
columns_2 = ['TARGETID', 'COADD_FIBERSTATUS', 'TARGET_RA', 'TARGET_DEC', 'MORPHTYPE', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'PARALLAX', 'EBV', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_Z', 'MASKBITS', 'PHOTSYS', 'DESI_TARGET', 'BGS_TARGET', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'COADD_NUMTILE', 'SCND_TARGET']
columns_4 = ['TARGETID', 'TSNR2_ELG', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG']
emline_columns = ['TARGETID', 'OII_FLUX', 'OII_FLUX_IVAR', 'OIII_FLUX', 'OIII_FLUX_IVAR', 'HALPHA_FLUX', 'HALPHA_FLUX_IVAR', 'HBETA_FLUX', 'HBETA_FLUX_IVAR', 'HGAMMA_FLUX', 'HGAMMA_FLUX_IVAR', 'HDELTA_FLUX', 'HDELTA_FLUX_IVAR']

cat_stack = []

fns = glob.glob('/global/cfs/cdirs/desi/users/rongpu/spectro/cosmos_lowz/daily/healpix/redrock*.fits')
print(len(fns))

for fn in fns:

    tmp1 = Table(fitsio.read(fn, ext=1, columns=columns_1))
    tmp2 = Table(fitsio.read(fn, ext=2, columns=columns_2))
    tmp4 = Table(fitsio.read(fn, ext=4, columns=columns_4))

    emline_fn = fn.replace('redrock-', 'emline-')
    tmp5 = Table(fitsio.read(emline_fn, ext=1, columns=(emline_columns)))

    if not (np.all(tmp1['TARGETID']==tmp2['TARGETID']) and np.all(tmp1['TARGETID']==tmp4['TARGETID']) and np.all(tmp1['TARGETID']==tmp5['TARGETID'])):
        raise ValueError
    cat = tmp1.copy()
    cat = join(cat, tmp2, keys='TARGETID')
    cat = join(cat, tmp4, keys='TARGETID')
    cat = join(cat, tmp5, keys='TARGETID')
    cat['fn'] = os.path.basename(fn)
    cat_stack.append(cat)

cat = vstack(cat_stack)
print(len(cat), len(np.unique(cat['TARGETID'])))

# Select science targets
mask = cat['SCND_TARGET'] & 2**62>0
print(np.sum(mask), np.sum(mask)/len(mask))
cat = cat[mask]

# # Remove FIBERSTATUS!=0 fibers
# mask = cat['COADD_FIBERSTATUS']==0
# print('FIBERSTATUS',np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
# cat = cat[mask]

# # Remove "no data" fibers
# mask = cat['ZWARN'] & 2**9==0
# print('No data', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
# cat = cat[mask]

cat.rename_column('TARGETID', 'TARGETID_TERTIARY')
mask = ~np.in1d(cat['TARGETID_TERTIARY'], parent['TARGETID_TERTIARY'])
print(np.sum(mask))  # sanity check

cat.remove_columns(['EBV', 'FIBERFLUX_Z', 'FLUX_G', 'FLUX_R', 'FLUX_W1', 'FLUX_W2', 'FLUX_Z', 'MASKBITS', 'MORPHTYPE'])
cat = join(cat, parent, keys='TARGETID_TERTIARY', join_type='left').filled(0)

cat.write('/global/cfs/cdirs/desi/users/rongpu/data/desi2/cosmos_lowz/spectro/daily/cosmos_lowz_healpix_coadds.fits', overwrite=True)
