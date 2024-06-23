# Select a gmag-limited ELG sample and a brighter ELG sample

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


def select_extended_lrg(cat):

    cat = cat.copy()

    w1shift = 0.21

    mask_quality = np.full(len(cat), True)

    mask_quality &= (cat['FLUX_IVAR_R'] > 0) & (cat['FLUX_R'] > 0)   # ADM quality in r.
    mask_quality &= (cat['FLUX_IVAR_Z'] > 0) & (cat['FLUX_Z'] > 0) & (cat['FIBERFLUX_Z'] > 0)   # ADM quality in z.
    mask_quality &= (cat['FLUX_IVAR_W1'] > 0) & (cat['FLUX_W1'] > 0)  # ADM quality in W1.

    mask_quality &= (cat['GAIA_PHOT_G_MEAN_MAG'] == 0) | (cat['GAIA_PHOT_G_MEAN_MAG'] > 18)  # remove bright GAIA sources

    # ADM remove stars with zfibertot < 17.5 that are missing from GAIA.
    mask_quality &= cat['FIBERTOTFLUX_Z'] < 10**(-0.4*(17.5-22.5))

    # ADM observed in every band.
    mask_quality &= (cat['NOBS_G'] > 0) & (cat['NOBS_R'] > 0) & (cat['NOBS_Z'] > 0)

    # Apply masks
    maskbits = [1, 12, 13]
    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    # print(np.sum(~mask_clean)/len(mask_clean))
    mask_quality &= mask_clean

    gmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_G']*10**(0.4*3.214*cat['EBV']), 1e-7, None))
    rmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_R']*10**(0.4*2.165*cat['EBV']), 1e-7, None))
    zmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_Z']*10**(0.4*1.211*cat['EBV']), 1e-7, None))
    w1mag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_W1']*10**(0.4*0.184*cat['EBV']), 1e-7, None))
    zfibermag = 22.5 - 2.5*np.log10(np.clip(cat['FIBERFLUX_Z']*10**(0.4*1.211*cat['EBV']), 1e-7, None))

    gr = gmag - rmag
    rz = rmag - zmag
    zw1 = zmag - w1mag
    rw1 = rmag - w1mag

    mask1 = mask_quality.copy()
    mask2 = mask_quality.copy()

    if field=='south':
        mask1 &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
        mask1 &= zfibermag < 21.6                   # faint limit
        mask1 &= (gmag - w1mag > 2.9) | (rmag - w1mag > 1.8)  # low-z cuts
        mask1 &= (
            ((rmag - w1mag > ((w1mag - w1shift) - 17.14) * 1.8)
             & (rmag - w1mag > ((w1mag - w1shift) - 16.33) * 1.))
            | (rmag - w1mag > 3.2)
        )  # double sliding cuts and high-z extension
    else:
        mask1 &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
        mask1 &= zfibermag < 21.61                   # faint limit
        mask1 &= (gmag - w1mag > 2.97) | (rmag - w1mag > 1.8)  # low-z cuts
        mask1 &= (
            ((rmag - w1mag > ((w1mag - w1shift) - 17.13) * 1.83)
             & (rmag - w1mag > ((w1mag - w1shift) - 16.31) * 1.))
            | (rmag - w1mag > 3.3)
        )  # double sliding cuts and high-z extension

    mask2 &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
    if field=='south':
        mask2 &= zfibermag < 21.6                  # faint limit
    else:
        mask2 &= zfibermag < 21.61                   # faint limit
    mask2 &= (((rz-(-0.2))**2 + (gr-(2.3))**2)>2**2) & (rz>0.4)

    mask_lrg = mask1 | mask2

    return mask_lrg


def get_sample(sweep_fn):

    print(sweep_fn)

    sweep_path = os.path.join(sweep_dir, sweep_fn)

    cat = Table(fitsio.read(sweep_path, columns=columns))
    cat['TARGETID'] = encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE'])
    extended_lrg = select_extended_lrg(cat)
    cat = cat[extended_lrg]

    return cat


cat_stack = []

for field in ['north', 'south']:

    sweep_dir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0'.format(field)

    sweep_all_path = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))
    sweep_fn_list = [os.path.basename(sweep_all_path[ii]) for ii in range(len(sweep_all_path))]

    columns = ['OBJID', 'BRICKID', 'RELEASE', 'RA', 'DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FIBERFLUX_Z', 'FIBERTOTFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG', 'EBV', 'MASKBITS', 'NOBS_G', 'NOBS_R', 'NOBS_Z']

    print('Start!')

    time_start = time.time()

    # start multiple worker processes
    n_processess = 128
    with Pool(processes=n_processess) as pool:
        res = pool.map(get_sample, sweep_fn_list)

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
cat_stack.write('/global/cfs/cdirs/desicollab/users/rongpu/tmp/lrg_desi_1b_sample_new.fits')

