# Select ELG targets with DESI-corrected EBV map

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp
from multiprocessing import Pool

from desitarget.targets import decode_targetid, encode_targetid

sys.path.append(os.path.expanduser('~/git/desi-targets/dr9/create_target_catalogs/main/'))
from select_desi_targets import select_elg


def get_elgs(sweep_fn):

    sweep_path = os.path.join(sweep_dir, sweep_fn)

    cat = Table(fitsio.read(sweep_path, columns=columns))
    cat['HPXPIXEL'] = hp.ang2pix(nside, cat['RA'], cat['DEC'], nest=False, lonlat=True)
    mask = np.in1d(cat['HPXPIXEL'], delta_gr_map['HPXPIXEL'])
    print(np.sum(mask)/len(mask))
    if np.sum(mask)==0:
        return None
    cat = cat[mask]
    cat = join(cat, delta_gr_map, keys='HPXPIXEL')
    cat['TARGETID'] = encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE'])

    cat['EBV_SFD'] = cat['EBV'].copy()
    cat['EBV'] = cat['EBV_SFD'] + (3.214-2.165) * cat['delta_gr_mean']
    cat['elglop'], cat['elgvlo'] = select_elg(cat)
    mask = cat['elglop'] | cat['elgvlo']
    print(np.sum(mask)/len(mask))
    if np.sum(mask)==0:
        return None
    cat = cat[mask]
    cat = cat[['TARGETID', 'RA', 'DEC', 'EBV', 'EBV_SFD', 'elglop', 'elgvlo']]

    return cat


nside = 64
delta_gr_map = Table(fitsio.read('/pscratch/sd/r/rongpu/ebv/desi_std/delta_gr_sv1sv3main_nside_{}.fits'.format(nside)))
print(len(delta_gr_map))
mask = delta_gr_map['n_star']>=3
delta_gr_map = delta_gr_map[mask]
print(len(delta_gr_map))
delta_gr_map = delta_gr_map[['HPXPIXEL', 'delta_gr_mean']]

cat_stack = []

for field in ['south', 'north']:

    sweep_dir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0'.format(field)

    sweep_all_path = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))
    sweep_fn_list = [os.path.basename(sweep_all_path[ii]) for ii in range(len(sweep_all_path))]

    columns = ['OBJID', 'BRICKID', 'RELEASE', 'RA', 'DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FIBERFLUX_G', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'EBV', 'MASKBITS', 'NOBS_G', 'NOBS_R', 'NOBS_Z']

    print('Start!')

    time_start = time.time()

    # start multiple worker processes
    n_processess = 128
    with Pool(processes=n_processess) as pool:
        res = pool.map(get_elgs, sweep_fn_list)

    # Remove None elements from the list
    for index in range(len(res)-1, -1, -1):
        if res[index] is None:
            res.pop(index)

    cat = vstack(res)

    if field=='north':
        mask_ns = (cat['DEC']>32.375)
    else:
        mask_ns = ((cat['DEC']<=32.375) | (cat['RA']<104) | (cat['RA']>280))
    cat = cat[mask_ns]
    cat['PHOTSYS'] = field[0].upper()
    print('Final combined catalog:', len(cat))

    cat_stack.append(cat)

cat_stack = vstack(cat_stack)
cat_stack.write('/pscratch/sd/r/rongpu/ebv/elg_targets_desi_ebv.fits', overwrite=True)

