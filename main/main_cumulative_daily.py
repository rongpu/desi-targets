############### The script is unfinished and does NOT work!!! ###############

# Get all Dark time LRGs, ELGs and QSOs
# python main_cumulative_all_dark.py LRG 20220226

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

# target_class: "LRG", "ELG", "QSO" or "BGS_ANY"
target_class = str(sys.argv[1])
target_class = target_class.upper()
date = str(sys.argv[2])

# The following target bits are the same in both main and SV3
target_bits = {'LRG': 0, 'ELG': 1, 'QSO': 2, 'BGS_ANY': 60}
target_bit = target_bits[target_class]

if target_class!='BGS_ANY':
    faflavor = 'maindark'
else:
    faflavor = 'mainbright'

tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
mask = tiles['FAFLAVOR']==faflavor
tiles = tiles[mask]
print(len(tiles), len(np.unique(tiles['TILEID'])))

mask = tiles['QA']=='good'
tiles = tiles[mask]
print(len(tiles), len(np.unique(tiles['TILEID'])))

columns_1 = ['TARGETID', 'CHI2', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'DELTACHI2']
columns_2 = ['TARGETID', 'PETAL_LOC', 'DEVICE_LOC', 'LOCATION', 'FIBER', 'COADD_FIBERSTATUS', 'TARGET_RA', 'TARGET_DEC', 'MORPHTYPE', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'PARALLAX', 'EBV', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'MASKBITS', 'PHOTSYS', 'DESI_TARGET', 'BGS_TARGET', 'TILEID', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'COADD_NUMTILE', 'FIBERASSIGN_X', 'FIBERASSIGN_Y']
columns_4 = ['TARGETID', 'TSNR2_ELG', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG']

# Use the archived (validated) catalogs
data_dir = '/global/cfs/cdirs/desi/spectro/redux/daily/tiles/archive'

cat_stack = []

for index in range(len(tiles)):

    tileid, archivedate, qanight = tiles['TILEID'][index], tiles['ARCHIVEDATE'][index], tiles['QANIGHT'][index]
    print(tileid)

    fn_list = sorted(glob.glob(os.path.join(data_dir, str(tileid), archivedate, 'redrock-0-{}-thru{}.fits'.format(tileid, qanight))))
    if len(fn_list)==0:
        fn_list = sorted(glob.glob(os.path.join(data_dir, str(tileid), '*/redrock-*.fits')))
    for fn in fn_list:
        tmp1 = Table(fitsio.read(fn, ext=1, columns=columns_1))
        tmp2 = Table(fitsio.read(fn, ext=2, columns=columns_2))
        tmp4 = Table(fitsio.read(fn, ext=4, columns=columns_4))
        if not (np.all(tmp1['TARGETID']==tmp2['TARGETID']) and np.all(tmp1['TARGETID']==tmp4['TARGETID'])):
            raise ValueError
        tmp = tmp1.copy()
        tmp = join(tmp, tmp2, keys='TARGETID')
        tmp = join(tmp, tmp4, keys='TARGETID')

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

cat.write('/global/cfs/cdirs/desi/users/rongpu/spectro/daily/main_cumulative_{}_{}.fits'.format(target_class, date), overwrite=True)
