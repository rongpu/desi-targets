from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp



tiles = Table.read('/global/cfs/cdirs/desi/spectro/redux/daily/tiles-daily.csv')
print(len(tiles))

mask = tiles['PROGRAM']=='dark1b'
tiles = tiles[mask]
print(len(tiles))

tiles['TILEID'].min(), tiles['TILEID'].max()

mask = tiles['EFFTIME_SPEC']>850
print(np.sum(~mask))
tiles = tiles[mask]
print(len(tiles))


def get_catalog(index):

    tileid = tiles['TILEID'][index]
    lastnight = tiles['LASTNIGHT'][index]
    fns = glob.glob(f'/global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/{tileid}/{lastnight}/redrock-*.fits')

    cat_stack = []

    for fn in fns:
        columns_1 = ['TARGETID', 'CHI2', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'DELTACHI2']
        columns_2 = ['TARGETID', 'PETAL_LOC', 'DEVICE_LOC', 'LOCATION', 'FIBER', 'COADD_FIBERSTATUS', 'TARGET_RA', 'TARGET_DEC', 'MORPHTYPE', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'PARALLAX', 'EBV', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_Z', 'MASKBITS', 'PHOTSYS', 'DESI_TARGET', 'BGS_TARGET', 'TILEID', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'COADD_NUMTILE', 'FIBERASSIGN_X', 'FIBERASSIGN_Y', 'PRIORITY', 'OBJTYPE']
        columns_4 = ['TARGETID', 'TSNR2_ELG', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG']
        columns_emline = ['TARGETID', 'OII_FLUX', 'OII_FLUX_IVAR']

        tmp1 = Table(fitsio.read(fn, ext=1, columns=columns_1))
        tmp2 = Table(fitsio.read(fn, ext=2, columns=columns_2))
        tmp4 = Table(fitsio.read(fn, ext=4, columns=columns_4))
        assert (np.all(tmp1['TARGETID']==tmp2['TARGETID']) and np.all(tmp1['TARGETID']==tmp4['TARGETID']))
        cat = tmp1.copy()
        cat = join(cat, tmp2, keys='TARGETID')
        cat = join(cat, tmp4, keys='TARGETID')
        cat['LASTNIGHT'] = np.array(os.path.basename(os.path.dirname(fn)), dtype=int)

        fn_emline = fn.replace('redrock-', 'emline-')
        emline = Table(fitsio.read(fn_emline, columns=columns_emline))
        cat = join(cat, emline, keys='TARGETID')

        mask = cat['DESI_TARGET'] & 2**3 > 0  # LGE
        cat = cat[mask]
        if len(cat)==0:
            print('No LGEs found;', fn)
            return None
        
        cat_stack.append(cat)
        
    cat = vstack(cat_stack)

    return cat


from multiprocessing import Pool
n_process = 128
with Pool(processes=n_process) as pool:
    res = pool.map(get_catalog, np.arange(len(tiles)))

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

cat = vstack(res)
print(len(cat))
cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/desi-ext/dark1b_lge_daily_20260211.fits')
