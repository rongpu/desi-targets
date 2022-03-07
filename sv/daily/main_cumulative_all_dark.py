# Get all Dark time LRGs, ELGs and QSOs

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits


coadd_type = 'cumulative'

tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
print(len(tiles), len(np.unique(tiles['TILEID'])))

mask = tiles['SURVEY']=='main'
mask &= tiles['FAPRGRM']=='dark'
tiles = tiles[mask]
print(len(tiles))

mask = tiles['OBSSTATUS']=='obsend'
# mask = tiles['QA']=='good'
tiles = tiles[mask]
print(len(tiles))

columns_1 = ['TARGETID', 'CHI2', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'DELTACHI2']
columns_2 = ['TARGETID', 'PETAL_LOC', 'DEVICE_LOC', 'LOCATION', 'FIBER', 'COADD_FIBERSTATUS', 'TARGET_RA', 'TARGET_DEC', 'MORPHTYPE', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'PARALLAX', 'EBV', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_Z', 'MASKBITS', 'PHOTSYS', 'DESI_TARGET', 'BGS_TARGET', 'TILEID', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'COADD_NUMTILE', 'FIBERASSIGN_X', 'FIBERASSIGN_Y']
columns_4 = ['TARGETID', 'TSNR2_ELG', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG']

data_dir_everest = '/global/cfs/cdirs/desi/spectro/redux/everest/tiles/{}'.format(coadd_type)
data_dir_daily = '/global/cfs/cdirs/desi/spectro/redux/daily/tiles/{}'.format(coadd_type)

cat_stack = []

for index, tileid in enumerate(tiles['TILEID']):

    print(tileid)

    fn_list = sorted(glob.glob(os.path.join(data_dir_everest, str(tileid), '*/redrock-*.fits')))
    if len(fn_list)==0:
        fn_list = sorted(glob.glob(os.path.join(data_dir_daily, str(tileid), '*/redrock-*.fits')))
    for fn in fn_list:
        tmp1 = Table(fitsio.read(fn, ext=1, columns=columns_1))
        tmp2 = Table(fitsio.read(fn, ext=2, columns=columns_2))
        tmp4 = Table(fitsio.read(fn, ext=4, columns=columns_4))
        if not (np.all(tmp1['TARGETID']==tmp2['TARGETID']) and np.all(tmp1['TARGETID']==tmp4['TARGETID'])):
            raise ValueError
        tmp = tmp1.copy()
        tmp = join(tmp, tmp2, keys='TARGETID')
        tmp = join(tmp, tmp4, keys='TARGETID')

        # Select LRG, ELG and QSO targets
        target_bits = [0, 1, 2]
        mask = np.full(len(tmp), False)
        for target_bit in target_bits:
            mask |= tmp['DESI_TARGET'] & 2**target_bit > 0
        tmp = tmp[mask]

        # night = fn[fn.find('-thru2021')+5:-5]
        # tmp['night'] = int(night)
        tmp['LASTNIGHT'] = tiles['LASTNIGHT'][index]
        tmp['QA'] = tiles['QA'][index]

        cat_stack.append(tmp)

cat = vstack(cat_stack)
print()
print(len(cat))

cat.write('/global/cfs/cdirs/desi/users/rongpu/spectro/daily/main_{}_all_dark_20220226.fits'.format(coadd_type), overwrite=True)
