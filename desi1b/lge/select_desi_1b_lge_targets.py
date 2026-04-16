# Select LGE sample for DESI-ext
# srun -N 1 -C cpu -c 256 -t 04:00:00 -q interactive python select_lrgs_for_desi_extension.py

# To get LRG mask:
# cd /global/u2/r/rongpu/git/desi-examples/bright_star_mask
# srun -N 1 -C cpu -c 256 -t 04:00:00 -q interactive python read_pixel_bitmask.py --tracer lrg --input /global/cfs/cdirs/desicollab/users/rongpu/data/desi-ext/desi_1b_lge_targets.fits --output /global/cfs/cdirs/desicollab/users/rongpu/data/desi-ext/desi_1b_lge_targets_lrgmask_v1.1.fits.gz


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


def select_extended_lrg_option_1(cat):

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

    mask_lrg = mask_quality.copy()

    if field=='south':
        mask_lrg &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
        mask_lrg &= zfibermag < 21.6                   # faint limit
        mask_lrg &= (gmag - w1mag > 2.9) | (rmag - w1mag > 1.8)  # low-z cuts
        mask_lrg &= (
            ((rmag - w1mag > ((w1mag - 0.275) - 17.14) * 1.8)
             & (rmag - w1mag > ((w1mag - 0.275) - 16.33) * 1.))
            | (rmag - w1mag > 3.3 - 0.1)
        )  # double sliding cuts and high-z extension
        mask_lrg &= ~(
            ((rmag - w1mag > (w1mag - 17.14) * 1.8)
             & (rmag - w1mag > (w1mag - 16.33) * 1.))
            | (rmag - w1mag > 3.3)
        )  # exclude DESI-1 LRGs
    else:
        mask_lrg &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
        mask_lrg &= zfibermag < 21.61                   # faint limit
        mask_lrg &= (gmag - w1mag > 2.97) | (rmag - w1mag > 1.8)  # low-z cuts
        mask_lrg &= (
            ((rmag - w1mag > ((w1mag - 0.275) - 17.13) * 1.83)
             & (rmag - w1mag > ((w1mag - 0.275) - 16.31) * 1.))
            | (rmag - w1mag > 3.4 - 0.1)
        )  # double sliding cuts and high-z extension
        mask_lrg &= ~(
            ((rmag - w1mag > (w1mag - 17.13) * 1.83)
             & (rmag - w1mag > (w1mag - 16.31) * 1.))
            | (rmag - w1mag > 3.4)
        )  # exclude DESI-1 LRGs

    return mask_lrg


def get_sample(sweep_fn):

    print(sweep_fn)

    sweep_path = os.path.join(sweep_dir, sweep_fn)

    cat = Table(fitsio.read(sweep_path, columns=columns))
    cat['TARGETID'] = encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE'])
    mask = select_extended_lrg_option_1(cat)
    cat = cat[mask]

    return cat


cat_stack = []

for field in ['north', 'south']:

    sweep_dir = '/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0'.format(field)

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

main = Table(fitsio.read('/dvs_ro/cfs/cdirs/desicollab/users/rongpu/targets/dr9.0/1.1.1/resolve/dr9_lrg_1.1.1_basic.fits', columns=['TARGETID']))
print(len(main))
cat_stack['main_lrg'] = np.in1d(cat_stack['TARGETID'], main['TARGETID'])

cat_stack.write('/global/cfs/cdirs/desicollab/users/rongpu/data/desi-ext/desi_1b_lge_targets.fits')

