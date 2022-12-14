# Create combined SV1 cumulative LRG catalog
# The same script also works for 1x_depth, 4x_depth and low_speed coadds

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
# coadd_type = '1x_depth'
# coadd_type = '4x_depth'
# coadd_type = 'low_speed'

tiles = Table(fitsio.read('/global/cfs/cdirs/desi/survey/observations/SV1/sv1-tiles.fits'))
print(len(tiles))
print(list(np.unique(tiles['TARGETS'])))

mask = tiles['TARGETS']=='QSO+LRG'
tiles = tiles[mask]
print(len(tiles))


columns_1 = ['TARGETID', 'CHI2', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'DELTACHI2']
columns_2 = ['TARGETID', 'PETAL_LOC', 'DEVICE_LOC', 'LOCATION', 'FIBER', 'COADD_FIBERSTATUS', 'TARGET_RA', 'TARGET_DEC', 'MORPHTYPE', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'PARALLAX', 'EBV', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_Z', 'MASKBITS', 'PHOTSYS', 'SV1_DESI_TARGET', 'SV1_BGS_TARGET', 'DESI_TARGET', 'BGS_TARGET', 'TILEID', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'COADD_NUMTILE']
columns_4 = ['TARGETID', 'TSNR2_ELG', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG']

tileid_list = tiles['TILEID']

if coadd_type!='low_speed':
    data_dir = '/global/cfs/cdirs/desi/spectro/redux/everest/tiles/{}'.format(coadd_type)
else:
    data_dir = '/global/cfs/cdirs/desi/users/rongpu/everest/low_speed'

cat_stack = []

for tileid in tileid_list:

    print(tileid)

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
        
        mask = tmp['SV1_DESI_TARGET'] & 2**0 > 0
        tmp = tmp[mask]

        if coadd_type in ['1x_depth', '4x_depth', 'low_speed']:
            str_tmp = os.path.join(data_dir, str(tileid)) + '/'
            subset = int(fn[len(str_tmp):][:fn[len(str_tmp):].find('/redrock')])
            tmp['subset'] = subset
        
        cat_stack.append(tmp)

cat = vstack(cat_stack)
print()
print(len(cat))

############################ Add LRG mask ############################

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

cat.write('/global/cfs/cdirs/desi/users/rongpu/spectro/everest/sv1_{}_lrg.fits'.format(coadd_type), overwrite=False)
