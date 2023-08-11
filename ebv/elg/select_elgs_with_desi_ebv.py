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


def select_elg(cat):
    '''
    columns = ['OBJID', 'BRICKID', 'RELEASE', 'RA', 'DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FIBERFLUX_G', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'EBV', 'MASKBITS', 'NOBS_G', 'NOBS_R', 'NOBS_Z']
    fn = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/sweep-000p000-010p005.fits'
    cat = Table(fitsio.read(fn, columns=columns))
    elg_original, mask_elgvlo = select_elg(cat)
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
    gmag = 22.5 - 2.5 * np.log10(np.clip(cat['FLUX_G']*10**(0.4*3.214*cat['EBV_DESI']), 1e-7, None))
    rmag = 22.5 - 2.5 * np.log10(np.clip(cat['FLUX_R']*10**(0.4*2.165*cat['EBV_DESI']), 1e-7, None))
    zmag = 22.5 - 2.5 * np.log10(np.clip(cat['FLUX_Z']*10**(0.4*1.211*cat['EBV_DESI']), 1e-7, None))
    gfibermag = 22.5 - 2.5 * np.log10(np.clip(cat['FIBERFLUX_G']*10**(0.4*3.214*cat['EBV_DESI']), 1e-7, None))

    elg_original = mask_quality.copy()
    elg_original &= gmag > 20                       # bright cut.
    elg_original &= rmag - zmag > 0.15                  # blue cut.
    elg_original &= gfibermag < 24.1  # faint cut.
    elg_original &= gmag - rmag < 0.5*(rmag - zmag) + 0.1  # remove stars, low-z galaxies.
    # ADM high-priority OII flux cut.
    elg_original &= gmag - rmag < -1.2*(rmag - zmag) + 1.3

    elg_gmag = mask_quality.copy()
    elg_gmag &= gmag > 20                       # bright cut.
    elg_gmag &= rmag - zmag > 0.15                  # blue cut.
    elg_gmag &= gmag < 23.5  # faint cut.
    elg_gmag &= gmag - rmag < 0.5*(rmag - zmag) + 0.1  # remove stars, low-z galaxies.
    # ADM high-priority OII flux cut.
    elg_gmag &= gmag - rmag < -1.2*(rmag - zmag) + 1.3

    elg_brighter = mask_quality.copy()
    elg_brighter &= gmag > 20                       # bright cut.
    elg_brighter &= rmag - zmag > 0.15                  # blue cut.
    elg_brighter &= gfibermag < 23.6  # faint cut.
    elg_brighter &= gmag - rmag < 0.5*(rmag - zmag) + 0.1  # remove stars, low-z galaxies.
    # ADM high-priority OII flux cut.
    elg_brighter &= gmag - rmag < -1.2*(rmag - zmag) + 1.3

    elg_gmag_brighter = mask_quality.copy()
    elg_gmag_brighter &= gmag > 20                       # bright cut.
    elg_gmag_brighter &= rmag - zmag > 0.15                  # blue cut.
    elg_gmag_brighter &= gmag < 23.0  # faint cut.
    elg_gmag_brighter &= gmag - rmag < 0.5*(rmag - zmag) + 0.1  # remove stars, low-z galaxies.
    # ADM high-priority OII flux cut.
    elg_gmag_brighter &= gmag - rmag < -1.2*(rmag - zmag) + 1.3

    return elg_original, elg_gmag, elg_brighter, elg_gmag_brighter


def get_sample(sweep_fn):

    print(sweep_fn)

    sweep_path = os.path.join(sweep_dir, sweep_fn)

    cat = Table(fitsio.read(sweep_path, columns=columns))
    cat['HPXPIXEL'] = hp.ang2pix(nside, cat['RA'], cat['DEC'], nest=False, lonlat=True)
    mask = np.in1d(cat['HPXPIXEL'], delta_gr_map['HPXPIXEL'])
    print(np.sum(mask)/len(mask))
    if np.sum(mask)==0:
        return None
    cat = cat[mask]
    cat.rename_column('EBV', 'EBV_SFD')
    cat = join(cat, delta_gr_map, keys='HPXPIXEL')
    cat['TARGETID'] = encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE'])

    elg_original, elg_gmag, elg_brighter, elg_gmag_brighter = select_elg(cat)
    cat['elg_original'] = elg_original.copy()
    cat['elg_gmag'] = elg_gmag.copy()
    cat['elg_brighter'] = elg_brighter.copy()
    cat['elg_gmag_brighter'] = elg_gmag_brighter.copy()
    mask = elg_original | elg_gmag | elg_brighter | elg_gmag_brighter
    if np.sum(mask)==0:
        return None
    cat = cat[mask]
    cat = cat[['TARGETID', 'RA', 'DEC', 'EBV_DESI', 'EBV_SFD', 'n_star', 'elg_original', 'elg_gmag', 'elg_brighter', 'elg_gmag_brighter']]

    return cat


nside = 128
delta_gr_map = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/v0/desi_std/maps/delta_gr_map_all_{}.fits'.format(nside)))
print(len(delta_gr_map))
mask = delta_gr_map['n_star']>=1
delta_gr_map = delta_gr_map[mask]
print(len(delta_gr_map))
delta_gr_map.rename_column('EBV', 'EBV_SFD')
delta_gr_map['EBV_DESI'] = delta_gr_map['delta_gr_hlmean']/1.049
delta_gr_map = delta_gr_map[['HPXPIXEL', 'EBV_DESI', 'n_star']]

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
cat_stack.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/v0/targets/alternative_elg_targets_desi_ebv.fits.gz', overwrite=True)

