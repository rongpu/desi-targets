# Create combined SV1 cumulative LRG catalog

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

tiles = Table(fitsio.read('/global/cfs/cdirs/desi/survey/observations/SV1/sv1-tiles.fits'))
print(len(tiles))
print(list(np.unique(tiles['TARGETS'])))

mask = tiles['TARGETS']=='QSO+LRG'
tiles = tiles[mask]
print(len(tiles))


columns_1 = ['TARGETID', 'CHI2', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'DELTACHI2']
columns_2 = ['TARGETID', 'PETAL_LOC', 'DEVICE_LOC', 'LOCATION', 'FIBER', 'COADD_FIBERSTATUS', 'TARGET_RA', 'TARGET_DEC', 'MORPHTYPE', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'PARALLAX', 'EBV', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'MASKBITS', 'PHOTSYS', 'SV1_DESI_TARGET', 'SV1_BGS_TARGET', 'DESI_TARGET', 'BGS_TARGET', 'TILEID', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'COADD_NUMTILE']
columns_4 = ['TARGETID', 'TSNR2_ELG', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG']
columns_emline = ['TARGETID', 'OII_FLUX', 'OII_FLUX_IVAR', 'OIII_FLUX', 'OIII_FLUX_IVAR', 'HALPHA_FLUX', 'HALPHA_FLUX_IVAR', 'HBETA_FLUX', 'HBETA_FLUX_IVAR', 'HGAMMA_FLUX', 'HGAMMA_FLUX_IVAR', 'HDELTA_FLUX', 'HDELTA_FLUX_IVAR']

tileid_list = tiles['TILEID']

top_data_dir = '/global/cfs/cdirs/desi/spectro/redux'
data_dir = os.path.join(top_data_dir, 'fuji/tiles/{}'.format(coadd_type))

cat_stack = []

for tileid in tileid_list:

    print(tileid)

    # tile_dir = os.path.join(data_dir, str(tileid))
    # lastnight = glob.glob(os.path.join(tile_dir, '*'))
    # if len(lastnight)==1:
    #     lastnight = int(os.path.basename(lastnight[0]))
    # elif len(lastnight)>1:
    #     raise ValueError('More than one lastnight: '+tile_dir)
    # elif len(lastnight)==0:
    #     print('Does not exist: '+tile_dir)
    #     continue

    fn_list = sorted(glob.glob(os.path.join(data_dir, str(tileid), '*/redrock-*.fits')))

    for fn in fn_list:
        tmp1 = Table(fitsio.read(fn, ext=1, columns=columns_1))
        tmp2 = Table(fitsio.read(fn, ext=2, columns=columns_2))
        tmp4 = Table(fitsio.read(fn, ext=4, columns=columns_4))
        emline_fn = fn.replace('redrock-', 'emline-')
        tmp5 = Table(fitsio.read(emline_fn, ext=1, columns=(columns_emline)))

        if not (np.all(tmp1['TARGETID']==tmp2['TARGETID']) and np.all(tmp1['TARGETID']==tmp4['TARGETID']) and np.all(tmp1['TARGETID']==tmp5['TARGETID'])):
            raise ValueError
        tmp = tmp1.copy()
        tmp = join(tmp, tmp2, keys='TARGETID')
        tmp = join(tmp, tmp4, keys='TARGETID')
        tmp = join(tmp, tmp5, keys='TARGETID')

        mask = tmp['SV1_DESI_TARGET'] & 2**0 > 0
        tmp = tmp[mask]

        if len(tmp)==0:
            print('No LRGs: ', fn)
            continue

        # tmp['LASTNIGHT'] = lastnight
        tmp['fn'] = fn[len(top_data_dir)+1:]

        if coadd_type in ['1x_depth', '4x_depth', 'lowspeed']:
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

main_dir = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.1.1/resolve'
lrg = Table(fitsio.read(os.path.join(main_dir, 'dr9_lrg_1.1.1_basic.fits'), columns=['TARGETID']))

mask = np.in1d(cat['TARGETID'], lrg['TARGETID'])
cat['main_lrg'] = False
cat['main_lrg'][mask] = True

#############################################################################

cat.write('/global/cfs/cdirs/desi/users/rongpu/spectro/fugu/sv1_{}_lrg.fits'.format(coadd_type), overwrite=True)
