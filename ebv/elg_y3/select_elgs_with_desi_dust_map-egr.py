# Use DESI reddening map to re-select ELGs
# Try two kinds of extinction corrections:
# 1. Simply replace EBV_SFD with EBV_GR (g-r-based EBV)
# 2. Separate reddening corrections for g-r and r-z colors, with EBV_GR and EBV_RZ, respectively; use EBV_GR for the extinction correction of the magnitudes
# This script does option 1.
# Three dust map resolutions are used: 128, 256 and 512.

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

sys.path.append(os.path.expanduser('~/git/desi-targets/useful/'))
from select_desi_targets import select_elg


def get_sample(sweep_fn):

    print(sweep_fn)

    sweep_path = os.path.join(sweep_dir, sweep_fn)

    cat = Table(fitsio.read(sweep_path, columns=columns))
    cat['HPXPIXEL'] = hp.ang2pix(nside, cat['RA'], cat['DEC'], nest=False, lonlat=True)
    mask = np.in1d(cat['HPXPIXEL'], dust_map['HPXPIXEL'])
    print(np.sum(mask)/len(mask))
    if np.sum(mask)==0:
        return None
    cat = cat[mask]
    cat = join(cat, dust_map, keys='HPXPIXEL', join_type='left')

    cat['EBV'] = cat['EBV_GR'].copy()

    cat['TARGETID'] = encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE'])

    mask_elglop, mask_elgvlo = select_elg(cat)
    cat['elglop'] = mask_elglop.copy()
    cat['elgvlo'] = mask_elgvlo.copy()
    mask = mask_elglop | mask_elgvlo
    if np.sum(mask)==0:
        return None
    cat = cat[mask]
    cat = cat[['TARGETID', 'RA', 'DEC', 'elglop', 'elgvlo']]

    return cat


for nside in [128, 256, 512]:

    dust_map = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars_y3/v0.2/final_maps/desi_dust_gr_{}.fits'.format(nside)))
    print(len(dust_map))
    dust_map = dust_map[['HPXPIXEL', 'EBV_GR']]

    cat_stack = []

    for field in ['north', 'south']:

        sweep_dir = '/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0'.format(field)

        sweep_all_path = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))
        sweep_fn_list = [os.path.basename(sweep_all_path[ii]) for ii in range(len(sweep_all_path))]

        columns = ['OBJID', 'BRICKID', 'RELEASE', 'RA', 'DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FIBERFLUX_G', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'MASKBITS', 'NOBS_G', 'NOBS_R', 'NOBS_Z']

        print('Start!')

        time_start = time.time()

        # start multiple worker processes
        n_processess = 128
        with Pool(processes=n_processess) as pool:
            res = pool.map(get_sample, sweep_fn_list, chunksize=1)

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
    cat_stack.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/misc/desi_targets_with_desi_ebv_y3/elg_targets-desi_egr-{}.fits'.format(nside), overwrite=False)

