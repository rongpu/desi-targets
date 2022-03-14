# Create combined SV3 perexp LRG catalog

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

tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-sv3.ecsv')
print(len(tiles))
mask = tiles['PROGRAM']=='DARK'
tiles = tiles[mask]
print(len(tiles))

mask = tiles['STATUS']=='done'
tiles = tiles[mask]
print(len(tiles))

columns_1 = ['TARGETID', 'CHI2', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'DELTACHI2']
columns_2 = ['TARGETID', 'PETAL_LOC', 'DEVICE_LOC', 'LOCATION', 'FIBER', 'COADD_FIBERSTATUS', 'TARGET_RA', 'TARGET_DEC', 'MORPHTYPE', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'PARALLAX', 'EBV', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'MASKBITS', 'PHOTSYS', 'SV3_DESI_TARGET', 'SV3_BGS_TARGET', 'DESI_TARGET', 'BGS_TARGET', 'TILEID', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'COADD_NUMTILE']
columns_4 = ['TARGETID', 'TSNR2_ELG', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG']

tileid_list = tiles['TILEID']

top_data_dir = '/global/cfs/cdirs/desi/spectro/redux'
data_dir = os.path.join(top_data_dir, 'fuji/tiles/{}'.format(coadd_type))

cat_stack = []

for tileid in tileid_list:

    print(tileid)

    tile_dir = os.path.join(data_dir, str(tileid))
    lastnight = glob.glob(os.path.join(tile_dir, '*'))
    if len(lastnight)==1:
        lastnight = int(os.path.basename(lastnight[0]))
    elif len(lastnight)>1:
        raise ValueError('More than one lastnight: '+tile_dir)
    elif len(lastnight)==0:
        raise ValueError('Does not exist: '+tile_dir)

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

        mask = tmp['SV3_DESI_TARGET'] & 2**0 > 0
        tmp = tmp[mask]

        if len(tmp)==0:
            print('No LRGs: ', fn)
            continue

        tmp['LASTNIGHT'] = lastnight
        tmp['fn'] = fn[len(top_data_dir)+1:]

        cat_stack.append(tmp)

cat = vstack(cat_stack)
print()
print(len(cat))

############################ Add LRG mask ############################

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

############################ Add main LRG flag ############################

main_dir = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve'

lrg = []
for field in ['north', 'south']:
    lrg.append(Table(fitsio.read(os.path.join(main_dir, 'dr9_lrg_{}_1.0.0_basic.fits'.format(field)), columns=['TARGETID'])))
lrg = vstack(lrg)

mask = np.in1d(cat['TARGETID'], lrg['TARGETID'])
cat['main_lrg'] = False
cat['main_lrg'][mask] = True

#############################################################################

cat.write('/global/cfs/cdirs/desi/users/rongpu/spectro/fugu/sv3_{}_lrg.fits'.format(coadd_type), overwrite=True)
