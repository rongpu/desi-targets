# Example:
# salloc -N 1 -C cpu -q interactive -t 4:00:00
# python create_per_tracer_catalogs-elg.py

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack
import fitsio

from multiprocessing import Pool

from desitarget.targets import encode_targetid

n_processes = 128

sweep_dir_north = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/north/sweep/9.0'
sweep_dir_south = '/global/cfs/cdirs/desicollab/users/rongpu/data/dr9/sweep_symlinks'

output_dir = '/global/cfs/cdirs/desicollab/users/rongpu/targets/dr9.0/zp_offset_corrected'

basic_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'PHOTSYS', 'TARGETID', 'ELG_LOP', 'ELG_VLO']

photom_columns = ['TYPE', 'EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z',
                  'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'FLUX_W1', 'FLUX_W2',
                  'FLUX_IVAR_W1', 'FLUX_IVAR_W2', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2',
                  'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG']

more_columns_1 = ['GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG', 'GAIA_ASTROMETRIC_EXCESS_NOISE', 'FITBITS',
                  'FRACFLUX_G', 'FRACFLUX_R', 'FRACFLUX_Z', 'FRACFLUX_W1', 'FRACFLUX_W2', 'FRACMASKED_G', 'FRACMASKED_R',
                  'FRACMASKED_Z', 'FRACIN_G', 'FRACIN_R', 'FRACIN_Z', 'FIBERTOTFLUX_G',
                  'SHAPE_R', 'SHAPE_R_IVAR', 'SHAPE_E1', 'SHAPE_E2', 'SERSIC', 'DCHISQ']

more_columns_2 = ['GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z',
                  'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'PSFDEPTH_W1', 'PSFDEPTH_W2',
                  'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z']

sweep_columns_all = basic_columns + photom_columns + more_columns_1 + more_columns_2
sweep_columns_all = list(set(sweep_columns_all))  # unique columns

sweep_columns_all.remove('TARGETID')
sweep_columns_all.remove('PHOTSYS')
sweep_columns_all.remove('ELG_LOP')
sweep_columns_all.remove('ELG_VLO')
sweep_columns_all += ['OBJID', 'BRICKID', 'RELEASE']


def get_elgs(sweep_path):

    print(os.path.basename(sweep_path))

    # sweep_path = os.path.join(sweep_dir, sweep_fn)

    cat = Table(fitsio.read(sweep_path, columns=sweep_columns_all))

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
    gmag = 22.5 - 2.5 * np.log10(np.clip(cat['FLUX_G']*10**(0.4*3.214*cat['EBV']), 1e-7, None))
    rmag = 22.5 - 2.5 * np.log10(np.clip(cat['FLUX_R']*10**(0.4*2.165*cat['EBV']), 1e-7, None))
    zmag = 22.5 - 2.5 * np.log10(np.clip(cat['FLUX_Z']*10**(0.4*1.211*cat['EBV']), 1e-7, None))
    gfibermag = 22.5 - 2.5 * np.log10(np.clip(cat['FIBERFLUX_G']*10**(0.4*3.214*cat['EBV']), 1e-7, None))

    mask_elglop = mask_quality.copy()

    mask_elglop &= gmag > 20                       # bright cut.
    mask_elglop &= rmag - zmag > 0.15                  # blue cut.
    mask_elglop &= gfibermag < 24.1  # faint cut.
    mask_elglop &= gmag - rmag < 0.5*(rmag - zmag) + 0.1  # remove stars, low-z galaxies.

    mask_elgvlo = mask_elglop.copy()

    # ADM low-priority OII flux cut.
    mask_elgvlo &= gmag - rmag < -1.2*(rmag - zmag) + 1.6
    mask_elgvlo &= gmag - rmag >= -1.2*(rmag - zmag) + 1.3

    # ADM high-priority OII flux cut.
    mask_elglop &= gmag - rmag < -1.2*(rmag - zmag) + 1.3

    cat['ELG_LOP'] = mask_elglop.copy()
    cat['ELG_VLO'] = mask_elgvlo.copy()
    cat = cat[mask_elglop | mask_elgvlo]

    return cat


cat = []

for field in ['north', 'south']:

    if field=='north':
        sweep_dir = sweep_dir_north
    else:
        sweep_dir = sweep_dir_south

    sweep_path_list = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))

    with Pool(processes=n_processes) as pool:
        res = pool.map(get_elgs, sweep_path_list)

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

# cat_basic = cat[basic_columns].copy()
# cat_photom = cat[photom_columns].copy()
# cat_more_1 = cat[more_columns_1].copy()
# cat_more_2 = cat[more_columns_2].copy()

# cat_basic_path = os.path.join(output_dir, 'dr9_elg_basic.fits')
# cat_photom_path = os.path.join(output_dir, 'dr9_elg_photom.fits')
# cat_more_1_path = os.path.join(output_dir, 'dr9_elg_more_1.fits')
# cat_more_2_path = os.path.join(output_dir, 'dr9_elg_more_2.fits')

# cat_basic.write(cat_basic_path)
# cat_photom.write(cat_photom_path)
# cat_more_1.write(cat_more_1_path)
# cat_more_2.write(cat_more_2_path)

cat_basic = cat[basic_columns].copy()
cat_basic_path = os.path.join(output_dir, 'dr9_elg_basic.fits')
cat_basic.write(cat_basic_path)
