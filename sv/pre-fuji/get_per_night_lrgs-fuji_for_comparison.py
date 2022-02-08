# Compare with first Fuji data; [desi-data 5835]

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits


coadd_type = 'pernight'

################################################## SV1 ##################################################

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/spectro/pre-fuji/f5/sv1_{}_lrg.fits'.format(coadd_type)))

columns_1 = ['TARGETID', 'CHI2', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'DELTACHI2']
columns_2 = ['TARGETID', 'PETAL_LOC', 'DEVICE_LOC', 'LOCATION', 'FIBER', 'COADD_FIBERSTATUS', 'TARGET_RA', 'TARGET_DEC', 'MORPHTYPE', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'PARALLAX', 'EBV', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_Z', 'MASKBITS', 'PHOTSYS', 'SV1_DESI_TARGET', 'SV1_BGS_TARGET', 'DESI_TARGET', 'BGS_TARGET', 'TILEID', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'COADD_NUMTILE']
columns_4 = ['TARGETID', 'TSNR2_ELG', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG']

data_dir = '/global/cfs/cdirs/desi/spectro/redux/fuji/tiles/{}'.format(coadd_type)

cat_stack = []

for tileid in np.unique(cat['TILEID']):

    print(tileid)

    mask = cat['TILEID']==tileid
    nights = np.unique(cat['NIGHT'][mask])

    for night in nights:

        fn_list = sorted(glob.glob(os.path.join(data_dir, str(tileid), str(night), 'redrock-*.fits')))
        for fn in fn_list:
            tmp1 = Table(fitsio.read(fn, ext=1, columns=columns_1))
            tmp2 = Table(fitsio.read(fn, ext=2, columns=columns_2))
            tmp4 = Table(fitsio.read(fn, ext=4, columns=columns_4))
            if not (np.all(tmp1['TARGETID']==tmp2['TARGETID']) and np.all(tmp1['TARGETID']==tmp4['TARGETID'])):
                raise ValueError
            tmp = tmp1.copy()
            tmp = join(tmp, tmp2, keys='TARGETID')
            tmp = join(tmp, tmp4, keys='TARGETID')

            mask = tmp['SV1_DESI_TARGET'] & 2**0 > 0
            tmp = tmp[mask]

            tmp['NIGHT'] = night

            ##############
            tmp['fn'] = fn
            ##############

            cat_stack.append(tmp)

if len(cat_stack)>0:

    cat = vstack(cat_stack)
    print()
    print(len(cat))

    # Add LRG mask
    sv1_dir = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/0.49.0'
    lrg = []
    for field in ['north', 'south']:
        tmp = Table(fitsio.read(os.path.join(sv1_dir, 'dr9_sv1_lrg_{}_0.49.0_basic.fits'.format(field)), columns=['TARGETID']))
        tmp1 = Table(fitsio.read(os.path.join(sv1_dir, 'dr9_sv1_lrg_{}_0.49.0_lrgmask_v1.fits'.format(field))))
        lrg.append(hstack([tmp, tmp1]))
    lrg = vstack(lrg)
    mask = np.in1d(lrg['TARGETID'], cat['TARGETID'])
    lrg = lrg[mask]
    cat = join(cat, lrg, keys='TARGETID')

    cat.write('/global/cfs/cdirs/desi/users/rongpu/spectro/pre-fuji/fuji/sv1_{}_lrg.fits'.format(coadd_type), overwrite=False)

################################################## SV3 ##################################################

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/spectro/pre-fuji/f5/sv3_{}_lrg.fits'.format(coadd_type)))

columns_1 = ['TARGETID', 'CHI2', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'DELTACHI2']
columns_2 = ['TARGETID', 'PETAL_LOC', 'DEVICE_LOC', 'LOCATION', 'FIBER', 'COADD_FIBERSTATUS', 'TARGET_RA', 'TARGET_DEC', 'MORPHTYPE', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'PARALLAX', 'EBV', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_Z', 'MASKBITS', 'PHOTSYS', 'SV3_DESI_TARGET', 'SV3_BGS_TARGET', 'DESI_TARGET', 'BGS_TARGET', 'TILEID', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'COADD_NUMTILE']
columns_4 = ['TARGETID', 'TSNR2_ELG', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG']

data_dir = '/global/cfs/cdirs/desi/spectro/redux/fuji/tiles/{}'.format(coadd_type)

cat_stack = []

for tileid in np.unique(cat['TILEID']):

    print(tileid)

    mask = cat['TILEID']==tileid
    nights = np.unique(cat['NIGHT'][mask])

    for night in nights:

        fn_list = sorted(glob.glob(os.path.join(data_dir, str(tileid), str(night), 'redrock-*.fits')))
        for fn in fn_list:
            tmp1 = Table(fitsio.read(fn, ext=1, columns=columns_1))
            tmp2 = Table(fitsio.read(fn, ext=2, columns=columns_2))
            tmp4 = Table(fitsio.read(fn, ext=4, columns=columns_4))
            if not (np.all(tmp1['TARGETID']==tmp2['TARGETID']) and np.all(tmp1['TARGETID']==tmp4['TARGETID'])):
                raise ValueError
            tmp = tmp1.copy()
            tmp = join(tmp, tmp2, keys='TARGETID')
            tmp = join(tmp, tmp4, keys='TARGETID')

            mask = tmp['SV3_DESI_TARGET'] & 2**0 > 0
            tmp = tmp[mask]

            tmp['NIGHT'] = night

            ##############
            tmp['fn'] = fn
            ##############

            cat_stack.append(tmp)

if len(cat_stack)>0:

    cat = vstack(cat_stack)
    print()
    print(len(cat))

    # Add LRG mask
    sv3_dir = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/0.57.0'
    lrg = []
    for field in ['north', 'south']:
        tmp = Table(fitsio.read(os.path.join(sv3_dir, 'dr9_sv3_lrg_{}_0.57.0_basic.fits'.format(field)), columns=['TARGETID']))
        tmp1 = Table(fitsio.read(os.path.join(sv3_dir, 'dr9_sv3_lrg_{}_0.57.0_lrgmask_v1.fits'.format(field))))
        lrg.append(hstack([tmp, tmp1]))
    lrg = vstack(lrg)
    mask = np.in1d(lrg['TARGETID'], cat['TARGETID'])
    lrg = lrg[mask]
    cat = join(cat, lrg, keys='TARGETID')

    cat.write('/global/cfs/cdirs/desi/users/rongpu/spectro/pre-fuji/fuji/sv3_{}_lrg.fits'.format(coadd_type), overwrite=False)

################################################## Main ##################################################

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/spectro/pre-fuji/f5/main_{}_lrg.fits'.format(coadd_type)))

columns_1 = ['TARGETID', 'CHI2', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'DELTACHI2']
columns_2 = ['TARGETID', 'PETAL_LOC', 'DEVICE_LOC', 'LOCATION', 'FIBER', 'COADD_FIBERSTATUS', 'TARGET_RA', 'TARGET_DEC', 'MORPHTYPE', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'PARALLAX', 'EBV', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_Z', 'MASKBITS', 'PHOTSYS', 'DESI_TARGET', 'BGS_TARGET', 'TILEID', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'COADD_NUMTILE', 'FIBERASSIGN_X', 'FIBERASSIGN_Y']
columns_4 = ['TARGETID', 'TSNR2_ELG', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG']

data_dir = '/global/cfs/cdirs/desi/spectro/redux/guadalupe/tiles/{}'.format(coadd_type)

cat_stack = []

for tileid in np.unique(cat['TILEID']):

    print(tileid)

    mask = cat['TILEID']==tileid
    nights = np.unique(cat['NIGHT'][mask])

    for night in nights:

        fn_list = sorted(glob.glob(os.path.join(data_dir, str(tileid), str(night), 'redrock-*.fits')))
        for fn in fn_list:
            tmp1 = Table(fitsio.read(fn, ext=1, columns=columns_1))
            tmp2 = Table(fitsio.read(fn, ext=2, columns=columns_2))
            tmp4 = Table(fitsio.read(fn, ext=4, columns=columns_4))
            if not (np.all(tmp1['TARGETID']==tmp2['TARGETID']) and np.all(tmp1['TARGETID']==tmp4['TARGETID'])):
                raise ValueError
            tmp = tmp1.copy()
            tmp = join(tmp, tmp2, keys='TARGETID')
            tmp = join(tmp, tmp4, keys='TARGETID')

            mask = tmp['DESI_TARGET'] & 2**0 > 0
            tmp = tmp[mask]

            tmp['NIGHT'] = night

            ##############
            tmp['fn'] = fn
            ##############

            cat_stack.append(tmp)

if len(cat_stack)>0:

    cat = vstack(cat_stack)
    print()
    print(len(cat))

    # Add LRG mask
    main_dir = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve'
    lrg = []
    for field in ['north', 'south']:
        tmp = Table(fitsio.read(os.path.join(main_dir, 'dr9_lrg_{}_1.0.0_basic.fits'.format(field)), columns=['TARGETID']))
        tmp1 = Table(fitsio.read(os.path.join(main_dir, 'dr9_lrg_{}_1.0.0_lrgmask_v1.fits'.format(field))))
        lrg.append(hstack([tmp, tmp1]))
    lrg = vstack(lrg)
    mask = np.in1d(lrg['TARGETID'], cat['TARGETID'])
    lrg = lrg[mask]
    cat = join(cat, lrg, keys='TARGETID')

    cat.write('/global/cfs/cdirs/desi/users/rongpu/spectro/pre-fuji/fuji/main_{}_lrg.fits'.format(coadd_type), overwrite=False)
