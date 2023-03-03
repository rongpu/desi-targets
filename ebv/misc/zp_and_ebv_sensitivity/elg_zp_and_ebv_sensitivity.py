# Get ELG targets with shifted ZP or EBV

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

    print(sweep_fn)

    sweep_path = os.path.join(sweep_dir, sweep_fn)

    cat = Table(fitsio.read(sweep_path, columns=columns))
    cat['TARGETID'] = encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE'])

    flags = Table()

    mask_all = np.full(len(cat), False)
    mask, _ = select_elg(cat)
    mask_all |= mask
    flags['original'] = mask.copy()

    for dmag in [-0.02, -0.01, 0.01, 0.02]:
        for band in ['g', 'r', 'z']:
            cat['FLUX_'+band.upper()] *= 10**(-0.4*(dmag))
            if band=='g':
                cat['FIBERFLUX_'+band.upper()] *= 10**(-0.4*(dmag))
            mask, _ = select_elg(cat)
            mask_all |= mask
            flags['dmag_{}_{}'.format(band, dmag)] = mask.copy()
            # restore values
            cat['FLUX_'+band.upper()] /= 10**(-0.4*(dmag))
            if band=='g':
                cat['FIBERFLUX_'+band.upper()] /= 10**(-0.4*(dmag))

    cat.rename_column('EBV', 'EBV_SFD')

    for ebv_additive in [-0.02, -0.01, -0.005, 0.005, 0.01, 0.02]:
        cat['EBV'] = cat['EBV_SFD'] + ebv_additive
        mask, _ = select_elg(cat)
        mask_all |= mask
        flags['ebv_add_{}'.format(ebv_additive)] = mask.copy()

    for ebv_multiplicative in [0.8, 0.9, 0.95, 1.05, 1.1, 1.2]:
        cat['EBV'] = cat['EBV_SFD'] * ebv_multiplicative
        mask, _ = select_elg(cat)
        mask_all |= mask
        flags['ebv_mult_{}'.format(ebv_multiplicative)] = mask.copy()

    if np.sum(mask_all)==0:
        return None

    cat = cat[mask_all]
    flags = flags[mask_all]

    cat = cat[['TARGETID', 'RA', 'DEC', 'EBV']]
    cat = hstack([cat, flags])

    return cat


cat_stack = []

for field in ['north', 'south']:

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
cat_stack.write('/pscratch/sd/r/rongpu/ebv/elg_lop_zp_and_ebv_sensitivity.fits')

