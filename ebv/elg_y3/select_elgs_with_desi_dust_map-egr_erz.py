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

# sys.path.append(os.path.expanduser('~/git/desi-targets/useful/'))
# from select_desi_targets import select_elg

r = {'g_south': 3.214, 'r_south': 2.165, 'z_south': 1.211, 'g_north': 3.258, 'r_north': 2.176, 'z_north': 1.199}

def select_elg(cat):
    '''
    columns = ['OBJID', 'BRICKID', 'RELEASE', 'RA', 'DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FIBERFLUX_G', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'EBV', 'MASKBITS', 'NOBS_G', 'NOBS_R', 'NOBS_Z']
    fn = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/sweep-000p000-010p005.fits'
    cat = Table(fitsio.read(fn, columns=columns))
    mask_elglop, mask_elgvlo = select_elg(cat)
    '''

    cat = cat.copy()
    cat.rename_columns(cat.colnames, [ii.upper() for ii in cat.colnames])

    mask_quality = np.full(len(cat), True)

    mask_quality &= (cat['FLUX_IVAR_G'] > 0) & (cat['FLUX_G'] > 0) & (cat['FIBERFLUX_G'] > 0)
    mask_quality &= (cat['FLUX_IVAR_R'] > 0) & (cat['FLUX_R'] > 0)
    mask_quality &= (cat['FLUX_IVAR_Z'] > 0) & (cat['FLUX_Z'] > 0)

    # ADM observed in every band.
    mask_quality &= (cat['NOBS_G'] > 0) & (cat['NOBS_R'] > 0) & (cat['NOBS_Z'] > 0)

    # Apply masks
    maskbits = [1, 12, 13]
    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    # print(np.sum(~mask_clean)/len(mask_clean))
    mask_quality &= mask_clean

    # gmag = 22.5 - 2.5 * np.log10((cat['FLUX_G'] / cat['MW_TRANSMISSION_G']).clip(1e-7))
    # rmag = 22.5 - 2.5 * np.log10((cat['FLUX_R'] / cat['MW_TRANSMISSION_R']).clip(1e-7))
    # zmag = 22.5 - 2.5 * np.log10((cat['FLUX_Z'] / cat['MW_TRANSMISSION_Z']).clip(1e-7))
    # gfibermag = 22.5 - 2.5 * np.log10((cat['FIBERFLUX_G'] / cat['MW_TRANSMISSION_G']).clip(1e-7))
    gmag = 22.5 - 2.5 * np.log10(np.clip(cat['FLUX_G']*10**(0.4*r['g_'+field]*cat['EBV']), 1e-7, None))
    rmag = 22.5 - 2.5 * np.log10(np.clip(cat['FLUX_R']*10**(0.4*r['r_'+field]*cat['EBV']), 1e-7, None))
    zmag = 22.5 - 2.5 * np.log10(np.clip(cat['FLUX_Z']*10**(0.4*r['z_'+field]*cat['EBV']), 1e-7, None))
    gfibermag = 22.5 - 2.5 * np.log10(np.clip(cat['FIBERFLUX_G']*10**(0.4*r['g_'+field]*cat['EBV']), 1e-7, None))

    gmag_gr = 22.5 - 2.5 * np.log10(np.clip(cat['FLUX_G']*10**(0.4*r['g_'+field]*cat['EBV_GR']), 1e-7, None))
    rmag_gr = 22.5 - 2.5 * np.log10(np.clip(cat['FLUX_R']*10**(0.4*r['r_'+field]*cat['EBV_GR']), 1e-7, None))
    rmag_rz = 22.5 - 2.5 * np.log10(np.clip(cat['FLUX_R']*10**(0.4*r['r_'+field]*cat['EBV_RZ']), 1e-7, None))
    zmag_rz = 22.5 - 2.5 * np.log10(np.clip(cat['FLUX_Z']*10**(0.4*r['z_'+field]*cat['EBV_RZ']), 1e-7, None))

    gr = gmag_gr - rmag_gr  # g-r color with g-r-based reddening correction
    rz = rmag_rz - zmag_rz  # r-z color with r-z-based reddening correction

    mask_elglop = mask_quality.copy()

    mask_elglop &= gmag > 20                       # bright cut.
    mask_elglop &= rz > 0.15                  # blue cut.
    mask_elglop &= gfibermag < 24.1  # faint cut.
    mask_elglop &= gr < 0.5*(rz) + 0.1  # remove stars, low-z galaxies.

    mask_elgvlo = mask_elglop.copy()

    # ADM low-priority OII flux cut.
    mask_elgvlo &= gr < -1.2*(rz) + 1.6
    mask_elgvlo &= gr >= -1.2*(rz) + 1.3

    # ADM high-priority OII flux cut.
    mask_elglop &= gr < -1.2*(rz) + 1.3

    return mask_elglop, mask_elgvlo


def get_sample(sweep_fn):

    print(sweep_fn)

    sweep_path = os.path.join(sweep_dir, sweep_fn)

    cat = Table(fitsio.read(sweep_path, columns=columns))
    cat['HPXPIXEL'] = hp.ang2pix(nside, cat['RA'], cat['DEC'], nest=False, lonlat=True)
    mask = np.in1d(cat['HPXPIXEL'], dust_map_gr['HPXPIXEL'])
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

    dust_map_gr = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars_y3/v0.2/final_maps/desi_dust_gr_{}.fits'.format(nside)))
    print(len(dust_map_gr))
    dust_map_gr = dust_map_gr[['HPXPIXEL', 'EBV_GR']]

    dust_map_rz = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars_y3/v0.2/final_maps/desi_dust_rz_{}.fits'.format(nside)))
    print(len(dust_map_rz))
    dust_map_rz = dust_map_rz[['HPXPIXEL', 'EBV_RZ']]

    # Combine the two maps
    idx = np.intersect1d(dust_map_gr['HPXPIXEL'], dust_map_rz['HPXPIXEL'])
    mask = np.in1d(dust_map_gr['HPXPIXEL'], idx)
    print(np.sum(mask), np.sum(mask)/len(mask))
    dust_map_gr = dust_map_gr[mask]
    mask = np.in1d(dust_map_rz['HPXPIXEL'], idx)
    print(np.sum(mask), np.sum(mask)/len(mask))
    dust_map_rz = dust_map_rz[mask]
    dust_map_gr.sort('HPXPIXEL')
    dust_map_rz.sort('HPXPIXEL')
    assert np.all(dust_map_gr['HPXPIXEL']==dust_map_rz['HPXPIXEL'])
    dust_map = dust_map_gr.copy()
    dust_map['EBV_RZ'] = dust_map_rz['EBV_RZ']

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
    cat_stack.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/misc/desi_targets_with_desi_ebv_y3/elg_targets-desi_egr_erz-{}.fits'.format(nside), overwrite=False)

