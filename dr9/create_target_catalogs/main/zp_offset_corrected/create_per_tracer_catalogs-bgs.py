from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack
import fitsio

from multiprocessing import Pool

from desitarget.targets import encode_targetid

n_processes = 128

sweep_dir_north = '/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/north/sweep/9.0'
sweep_dir_south = '/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/dr9/sweep_symlinks'

output_dir = '/global/cfs/cdirs/desicollab/users/rongpu/targets/dr9.0/zp_offset_corrected'

basic_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'PHOTSYS', 'TARGETID']

photom_columns = ['TYPE', 'EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z',
                  'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'FLUX_W1', 'FLUX_W2',
                  'FLUX_IVAR_W1', 'FLUX_IVAR_W2', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2',
                  'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG']

more_columns_1 = ['GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG', 'GAIA_ASTROMETRIC_EXCESS_NOISE', 'FITBITS',
                  'FRACFLUX_G', 'FRACFLUX_R', 'FRACFLUX_Z', 'FRACFLUX_W1', 'FRACFLUX_W2', 'FRACMASKED_G', 'FRACMASKED_R',
                  'FRACMASKED_Z', 'FRACIN_G', 'FRACIN_R', 'FRACIN_Z', 'FIBERTOTFLUX_G',
                  'SHAPE_R', 'SHAPE_R_IVAR', 'SHAPE_E1', 'SHAPE_E2', 'SERSIC', 'DCHISQ', 'REF_CAT',
                  'REF_EPOCH', 'PMRA', 'PMDEC']

more_columns_2 = ['GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z',
                  'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'PSFDEPTH_W1', 'PSFDEPTH_W2',
                  'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z']

sweep_columns_all = basic_columns + photom_columns + more_columns_1 + more_columns_2
sweep_columns_all = list(set(sweep_columns_all))  # unique columns

sweep_columns_all.remove('TARGETID')
sweep_columns_all.remove('PHOTSYS')
sweep_columns_all += ['OBJID', 'BRICKID', 'RELEASE']


def get_bgs(sweep_path):

    # sweep_path = os.path.join(sweep_dir, sweep_fn)

    cat = Table(fitsio.read(sweep_path, columns=sweep_columns_all))

    # Use undereddened Gaia and r-band fluxes
    grr = cat['GAIA_PHOT_G_MEAN_MAG'] - 22.5 + 2.5*np.log10(1e-16)
    ii = cat['FLUX_R'] > 0
    grr[ii] = cat['GAIA_PHOT_G_MEAN_MAG'][ii] - 22.5 + 2.5*np.log10(cat['FLUX_R'][ii])

    # Dereddening the fluxes
    cat['FLUX_G_EC'] = cat['FLUX_G'] * 10**(0.4*3.214*cat['EBV'])
    cat['FLUX_R_EC'] = cat['FLUX_R'] * 10**(0.4*2.165*cat['EBV'])
    cat['FLUX_Z_EC'] = cat['FLUX_Z'] * 10**(0.4*1.211*cat['EBV'])
    cat['FLUX_W1_EC'] = cat['FLUX_W1'] * 10**(0.4*0.184*cat['EBV'])
    cat['FIBERFLUX_R_EC'] = cat['FIBERFLUX_R'] * 10**(0.4*2.165*cat['EBV'])

    g = 22.5 - 2.5*np.log10(cat['FLUX_G_EC'].clip(1e-16))
    r = 22.5 - 2.5*np.log10(cat['FLUX_R_EC'].clip(1e-16))
    z = 22.5 - 2.5*np.log10(cat['FLUX_Z_EC'].clip(1e-16))
    w1 = 22.5 - 2.5*np.log10(cat['FLUX_W1_EC'].clip(1e-16))
    rfib = 22.5 - 2.5*np.log10(cat['FIBERFLUX_R_EC'].clip(1e-16))

    mask_quality = np.full(len(cat), True)

    mask_quality &= (cat['NOBS_G'] > 0) & (cat['NOBS_R'] > 0) & (cat['NOBS_Z'] > 0)
    mask_quality &= (cat['FLUX_IVAR_G'] > 0) & (cat['FLUX_IVAR_R'] > 0) & (cat['FLUX_IVAR_Z'] > 0)
    mask_quality &= ((grr > 0.6) | (cat['GAIA_PHOT_G_MEAN_MAG']==0))

    fmc = np.full(len(cat), False)
    fmc |= ((rfib < (2.9 + 1.2 + 1.0) + r) & (r < 17.8))
    fmc |= ((rfib < 22.9) & (r < 20.0) & (r > 17.8))
    fmc |= ((rfib < 2.9 + r) & (r > 20))
    mask_quality &= fmc

    # the SGA galaxies.
    mask_sga = np.array([(rc[0] == "L") if len(rc) > 0 else False for rc in cat['REF_CAT']])
    mask_quality |= mask_sga

    # Apply masks
    maskbits = [1, 13]
    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    # print(np.sum(~mask_clean)/len(mask_clean))
    mask_quality &= mask_clean

    mask_bgs = mask_quality.copy()

    if field=='south':
        mask_bgs &= cat['FLUX_R_EC'] > cat['FLUX_G_EC'] * 10**(-1.0/2.5)
        mask_bgs &= cat['FLUX_R_EC'] < cat['FLUX_G_EC'] * 10**(4.0/2.5)
        mask_bgs &= cat['FLUX_Z_EC'] > cat['FLUX_R_EC'] * 10**(-1.0/2.5)
        mask_bgs &= cat['FLUX_Z_EC'] < cat['FLUX_R_EC'] * 10**(4.0/2.5)
    else:
        mask_bgs &= cat['FLUX_R_EC'] > cat['FLUX_G_EC'] * 10**(-1.0/2.5)
        mask_bgs &= cat['FLUX_R_EC'] < cat['FLUX_G_EC'] * 10**(4.0/2.5)
        mask_bgs &= cat['FLUX_Z_EC'] > cat['FLUX_R_EC'] * 10**(-1.0/2.5)
        mask_bgs &= cat['FLUX_Z_EC'] < cat['FLUX_R_EC'] * 10**(4.0/2.5)

    # BASS r-mag offset with DECaLS.
    offset = 0.04

    mask_bgs_bright = mask_bgs.copy()
    if field=='south':
        mask_bgs_bright &= cat['FLUX_R_EC'] > 10**((22.5-19.5)/2.5)
        mask_bgs_bright &= cat['FLUX_R_EC'] <= 10**((22.5-12.0)/2.5)
        mask_bgs_bright &= cat['FIBERTOTFLUX_R'] <= 10**((22.5-15.0)/2.5)
    else:
        mask_bgs_bright &= cat['FLUX_R_EC'] > 10**((22.5-(19.5+offset))/2.5)
        mask_bgs_bright &= cat['FLUX_R_EC'] <= 10**((22.5-12.0)/2.5)
        mask_bgs_bright &= cat['FIBERTOTFLUX_R'] <= 10**((22.5-15.0)/2.5)

    mask_bgs_faint = mask_bgs.copy()
    if field=='south':
        mask_bgs_faint &= cat['FLUX_R_EC'] > 10**((22.5-20.175)/2.5)
        mask_bgs_faint &= cat['FLUX_R_EC'] <= 10**((22.5-19.5)/2.5)
        schlegel_color = (z - w1) - 3/2.5 * (g - r) + 1.2
        rfibcol = (rfib < 20.75) | ((rfib < 21.5) & (schlegel_color > 0.))
        mask_bgs_faint &= (rfibcol)
    else:
        mask_bgs_faint &= cat['FLUX_R_EC'] > 10**((22.5-(20.220))/2.5)
        mask_bgs_faint &= cat['FLUX_R_EC'] <= 10**((22.5-(19.5+offset))/2.5)
        schlegel_color = (z - w1) - 3/2.5 * (g - (r-offset)) + 1.2
        rfibcol = (rfib < 20.75+offset) | ((rfib < 21.5+offset) & (schlegel_color > 0.))
        mask_bgs_faint &= (rfibcol)

    cat['BGS_BRIGHT'] = mask_bgs_bright
    cat['BGS_FAINT'] = mask_bgs_faint
    mask = mask_bgs_bright | mask_bgs_faint
    cat = cat[mask]

    cat.remove_columns(['FLUX_G_EC', 'FLUX_R_EC', 'FLUX_Z_EC', 'FLUX_W1_EC', 'FIBERFLUX_R_EC'])

    return cat


cat = []

for field in ['north', 'south']:

    if field=='north':
        sweep_dir = sweep_dir_north
    else:
        sweep_dir = sweep_dir_south

    sweep_path_list = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))

    with Pool(processes=n_processes) as pool:
        res = pool.map(get_bgs, sweep_path_list)

    # Remove None elements from the list
    for index in range(len(res)-1, -1, -1):
        if res[index] is None:
            res.pop(index)

    tmp = vstack(res)

    if field=='north':
        tmp['PHOTSYS'] = 'N'
    else:
        tmp['PHOTSYS'] = 'S'

    cat.append(tmp)

cat = vstack(cat)

mask_north = (cat['PHOTSYS']=='N') & (cat['DEC']>32.375)
mask_south = (cat['PHOTSYS']=='S') & ((cat['DEC']<=32.375) | (cat['RA']<104) | (cat['RA']>280))
mask = mask_north | mask_south
cat = cat[mask]

cat['TARGETID'] = encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE'])

cat_basic = cat[basic_columns].copy()
cat_photom = cat[photom_columns].copy()
cat_more_1 = cat[more_columns_1].copy()
cat_more_2 = cat[more_columns_2].copy()

cat_basic_path = os.path.join(output_dir, 'dr9_bgs_basic.fits')
cat_photom_path = os.path.join(output_dir, 'dr9_bgs_photom.fits')
cat_more_1_path = os.path.join(output_dir, 'dr9_bgs_more_1.fits')
cat_more_2_path = os.path.join(output_dir, 'dr9_bgs_more_2.fits')

cat_basic.write(cat_basic_path)
cat_photom.write(cat_photom_path)
cat_more_1.write(cat_more_1_path)
cat_more_2.write(cat_more_2_path)

