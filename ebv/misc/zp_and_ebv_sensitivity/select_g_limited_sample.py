# Select g-band magnitude/fibermagnitude-limited samples

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


def g_selection(cat):
    '''
    columns = ['OBJID', 'BRICKID', 'RELEASE', 'RA', 'DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FIBERFLUX_G', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'EBV', 'MASKBITS', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'FIBERTOTFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG']
    fn = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/sweep-000p000-010p005.fits'
    cat = Table(fitsio.read(fn, columns=columns))
    mask_gmag, mask_gfiber = g_selection(cat)
    '''

    cat = cat.copy()
    cat.rename_columns(cat.colnames, [ii.upper() for ii in cat.colnames])

    mask_quality = np.full(len(cat), True)

    mask_quality &= (cat['FLUX_IVAR_G'] > 0) & (cat['FLUX_G'] > 0) & (cat['FIBERFLUX_G'] > 0)
    mask_quality &= (cat['FLUX_IVAR_R'] > 0) & (cat['FLUX_R'] > 0)
    mask_quality &= (cat['FLUX_IVAR_Z'] > 0) & (cat['FLUX_Z'] > 0)
    mask_quality &= (cat['FLUX_IVAR_W1'] > 0) & (cat['FLUX_W1'] > 0)  # ADM quality in W1.
    mask_quality &= (cat['GAIA_PHOT_G_MEAN_MAG'] == 0) | (cat['GAIA_PHOT_G_MEAN_MAG'] > 18)  # remove bright GAIA sources

    # ADM observed in every band.
    mask_quality &= (cat['NOBS_G'] > 0) & (cat['NOBS_R'] > 0) & (cat['NOBS_Z'] > 0)

    # ADM remove stars with zfibertot < 17.5 that are missing from GAIA.
    mask_quality &= cat['FIBERTOTFLUX_Z'] < 10**(-0.4*(17.5-22.5))

    # Apply masks
    maskbits = [1, 8, 9, 11, 12, 13]
    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    # print(np.sum(~mask_clean)/len(mask_clean))
    mask_quality &= mask_clean

    # gmag = 22.5 - 2.5 * np.log10((cat['FLUX_G'] / cat['MW_TRANSMISSION_G']).clip(1e-7))
    # rmag = 22.5 - 2.5 * np.log10((cat['FLUX_R'] / cat['MW_TRANSMISSION_R']).clip(1e-7))
    # zmag = 22.5 - 2.5 * np.log10((cat['FLUX_Z'] / cat['MW_TRANSMISSION_Z']).clip(1e-7))
    # gfibermag = 22.5 - 2.5 * np.log10((cat['FIBERFLUX_G'] / cat['MW_TRANSMISSION_G']).clip(1e-7))
    gmag = 22.5 - 2.5 * np.log10(np.clip(cat['FLUX_G']*10**(0.4*3.214*cat['EBV']), 1e-7, None))
    rmag = 22.5 - 2.5 * np.log10(np.clip(cat['FLUX_R']*10**(0.4*2.165*cat['EBV']), 1e-7, None))
    zmag = 22.5 - 2.5 * np.log10(np.clip(cat['FLUX_Z']*10**(0.4*1.211*cat['EBV']), 1e-7, None))
    w1mag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_W1']*10**(0.4*0.184*cat['EBV']), 1e-7, None))
    gfibermag = 22.5 - 2.5 * np.log10(np.clip(cat['FIBERFLUX_G']*10**(0.4*3.214*cat['EBV']), 1e-7, None))

    mask_quality &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.8  # non-stellar cut

    mask_gmag = mask_quality.copy()
    mask_gfiber = mask_quality.copy()

    mask_gmag &= gmag < 21.7  # faint cut.
    mask_gfiber &= gfibermag < 22.6  # faint cut.

    return mask_gmag, mask_gfiber


def get_sample(sweep_fn):

    print(sweep_fn)

    sweep_path = os.path.join(sweep_dir, sweep_fn)

    cat = Table(fitsio.read(sweep_path, columns=columns))
    cat['TARGETID'] = encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE'])
    mask_gmag, mask_gfiber = g_selection(cat)
    cat['gmag_sel'] = mask_gmag.copy()
    cat['gfiber_sel'] = mask_gfiber.copy()
    mask = mask_gmag | mask_gfiber
    cat = cat[mask]
    cat = cat[['TARGETID', 'RA', 'DEC', 'gmag_sel', 'gfiber_sel']]

    return cat


cat_stack = []

for field in ['north', 'south']:

    sweep_dir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0'.format(field)

    sweep_all_path = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))
    sweep_fn_list = [os.path.basename(sweep_all_path[ii]) for ii in range(len(sweep_all_path))]

    columns = ['OBJID', 'BRICKID', 'RELEASE', 'RA', 'DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FIBERFLUX_G', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'EBV', 'MASKBITS', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'FIBERTOTFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG']

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
cat_stack.write('/pscratch/sd/r/rongpu/ebv/g_limited_sample.fits')

