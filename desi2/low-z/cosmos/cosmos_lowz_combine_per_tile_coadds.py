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
columns_2 = ['TARGETID', 'PETAL_LOC', 'DEVICE_LOC', 'LOCATION', 'FIBER', 'COADD_FIBERSTATUS', 'TARGET_RA', 'TARGET_DEC', 'MORPHTYPE', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'PARALLAX', 'EBV', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_Z', 'MASKBITS', 'PHOTSYS', 'DESI_TARGET', 'BGS_TARGET', 'TILEID', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'COADD_NUMTILE', 'FIBERASSIGN_X', 'FIBERASSIGN_Y', 'SCND_TARGET', 'TILEID']
columns_4 = ['TARGETID', 'TSNR2_ELG', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG']

cat_stack = []
for tileid in [82637, 82638, 82639, 82640, 82641, 82642, 82643, 82644, 82645, 82646, 82647]:
    for petal in range(10):
        fn = glob.glob('/global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/{}/202204*/redrock-{}-{}-thru202204*.fits'.format(tileid, petal, tileid))[0]

        tmp1 = Table(fitsio.read(fn, ext=1, columns=columns_1))
        tmp2 = Table(fitsio.read(fn, ext=2, columns=columns_2))
        tmp4 = Table(fitsio.read(fn, ext=4, columns=columns_4))
        if not (np.all(tmp1['TARGETID']==tmp2['TARGETID']) and np.all(tmp1['TARGETID']==tmp4['TARGETID'])):
            raise ValueError
        cat = tmp1.copy()
        cat = join(cat, tmp2, keys='TARGETID')
        cat = join(cat, tmp4, keys='TARGETID')
        cat['fn'] = fn[60:]
        cat_stack.append(cat)

cat = vstack(cat_stack)
print(len(cat))

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

# # Remove "bad" fibers
# bad_fibers = [466, 500, 542, 552, 651, 817, 961, 1008, 1098, 1400, 1597, 2252, 2253, 2255, 2260, 2262, 2316, 2575, 2628, 2636, 2654, 2663, 2666, 2675, 2676, 2678, 2679, 2680, 2681, 2684, 2685, 2688, 2689, 2773, 3124, 3429, 3448, 3476, 3481, 3500, 3518, 3618, 3849, 3974, 3994, 4003, 4019, 4089, 4119, 4349, 4621, 4624, 4638, 4704, 4720, 4788, 4957, 4977]
# mask = ~np.in1d(cat['FIBER'], bad_fibers)
# print('Bad fibers', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
# cat = cat[mask]

cat.rename_column('TARGETID', 'TARGETID_TERTIARY')
mask = ~np.in1d(cat['TARGETID_TERTIARY'], parent['TARGETID_TERTIARY'])
print(np.sum(mask))  # sanity check

cat.remove_columns(['EBV', 'FIBERFLUX_Z', 'FLUX_G', 'FLUX_R', 'FLUX_W1', 'FLUX_W2', 'FLUX_Z', 'MASKBITS', 'MORPHTYPE'])
cat = join(cat, parent, keys='TARGETID_TERTIARY', join_type='left').filled(0)

cat.write('/global/cfs/cdirs/desi/users/rongpu/data/desi2/cosmos_lowz/spectro/daily/healpix/cosmos_lowz_per_tile_coadds.fits', overwrite=True)
