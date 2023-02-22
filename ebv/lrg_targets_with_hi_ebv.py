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
from select_desi_targets import select_lrg


def get_lrgs(sweep_fn):

    sweep_path = os.path.join(sweep_dir, sweep_fn)

    cat = Table(fitsio.read(sweep_path, columns=columns))
    # c = SkyCoord(cat['RA'], cat['DEC'], unit='deg').galactic
    # l, b = c.l.to_value('deg'), c.b.to_value('deg')
    # cat['pix'] = hp.ang2pix(nside, l, b, nest=False, lonlat=True)
    cat['pix'] = hp.ang2pix(nside, cat['RA'], cat['DEC'], nest=False, lonlat=True)
    mask = np.in1d(cat['pix'], ebv['pix'])
    print(np.sum(mask)/len(mask))
    if np.sum(mask)==0:
        return None
    cat = cat[mask]
    cat = join(cat, ebv, keys='pix')
    cat['TARGETID'] = encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE'])

    cat1 = cat.copy()
    cat1['EBV'] = cat['EBV_HI']
    mask = select_lrg(cat1, field=field)
    print(np.sum(mask)/len(mask))
    if np.sum(mask)==0:
        return None
    cat = cat[mask]
    cat = cat[['TARGETID', 'RA', 'DEC', 'EBV', 'EBV_HI']]

    return cat


cat_stack = []

for field in ['south', 'north']:

    sweep_dir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0'.format(field)

    sweep_all_path = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))
    sweep_fn_list = [os.path.basename(sweep_all_path[ii]) for ii in range(len(sweep_all_path))]

    # fn_hi = '/global/cfs/cdirs/desi/users/rongpu/useful/lenz_hi/ebv_lhd.hpx.fits'
    fn_hi = '/global/cfs/cdirs/desi/users/rongpu/useful/lenz_hi/ebv_lhd_equatorial.hpx.fits'
    ebv = Table(fitsio.read(fn_hi))
    nside = 1024
    ebv['pix'] = np.arange(hp.nside2npix(nside))
    ebv.rename_column('EBV', 'EBV_HI')

    mask = np.isfinite(ebv['EBV_HI'])
    ebv = ebv[mask]

    ebv['EBV_HI'] = ebv['EBV_HI'] / 0.884  # rescale to match the original SFD EBV

    columns = ['OBJID', 'BRICKID', 'RELEASE', 'RA', 'DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FIBERFLUX_Z', 'FIBERTOTFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG', 'EBV', 'MASKBITS', 'NOBS_G', 'NOBS_R', 'NOBS_Z']

    print('Start!')

    time_start = time.time()

    # start multiple worker processes
    n_processess = 128
    with Pool(processes=n_processess) as pool:
        res = pool.map(get_lrgs, sweep_fn_list)

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
cat_stack.write('/pscratch/sd/r/rongpu/ebv/lrg_targets_hi_ebv.fits')

