# Create the target catalog for main LRGs
# The resulting catalog should be the same as Adam's official LRG catalog

# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python create_lrg_catalog.py

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio
from astropy.io import fits
from multiprocessing import Pool

from desitarget.targets import encode_targetid


print('Start!')

time_start = time.time()


n_processes = 64

sweep_columns = ['RELEASE', 'BRICKID', 'OBJID', 'TYPE', 'RA', 'DEC', 'EBV',
'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2',
'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2',
'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2',
'FIBERFLUX_Z', 'FIBERTOTFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG',
'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'WISEMASK_W1']


def get_lrgs(sweep_path):

    # sweep_path = os.path.join(sweep_dir, sweep_fn)

    # cat = Table(fitsio.read(sweep_path, columns=sweep_columns))
    
    # Using astropy.io.fits solves the problem/bug in NERSC that makes it super slow to
    # read large FITS tables with fitsio in the interactive node
    hdu = fits.open(sweep_path)
    cat = Table()
    for col in sweep_columns:
        cat[col] = np.copy(hdu[1].data[col])

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
    print(np.sum(~mask_clean)/len(mask_clean))
    mask_quality &= mask_clean

    ######################### Selection cuts #########################
    gmag = 22.5 - 2.5 * np.log10((cat['FLUX_G'] / cat['MW_TRANSMISSION_G']).clip(1e-7))
    # ADM safe as these fluxes are set to > 0 in notinLRG_mask.
    rmag = 22.5 - 2.5 * np.log10((cat['FLUX_R'] / cat['MW_TRANSMISSION_R']).clip(1e-7))
    zmag = 22.5 - 2.5 * np.log10((cat['FLUX_Z'] / cat['MW_TRANSMISSION_Z']).clip(1e-7))
    w1mag = 22.5 - 2.5 * np.log10((cat['FLUX_W1'] / cat['MW_TRANSMISSION_W1']).clip(1e-7))
    zfibermag = 22.5 - 2.5 * np.log10((cat['FIBERFLUX_Z'] / cat['MW_TRANSMISSION_Z']).clip(1e-7))

    mask_lrg = np.full(len(cat), True)

    # if field=='south':
    if '/dr9/south/sweep/' in sweep_path:
        mask_lrg &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
        mask_lrg &= zfibermag < 21.6                   # faint limit
        mask_lrg &= (gmag - w1mag > 2.9) | (rmag - w1mag > 1.8)  # low-z cuts
        mask_lrg &= (
            ((rmag - w1mag > (w1mag - 17.14) * 1.8)
             & (rmag - w1mag > (w1mag - 16.33) * 1.))
            | (rmag - w1mag > 3.3)
        )  # double sliding cuts and high-z extension
    elif '/dr9/north/sweep/' in sweep_path:
        mask_lrg &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
        mask_lrg &= zfibermag < 21.61                   # faint limit
        mask_lrg &= (gmag - w1mag > 2.97) | (rmag - w1mag > 1.8)  # low-z cuts
        mask_lrg &= (
            ((rmag - w1mag > (w1mag - 17.13) * 1.83)
             & (rmag - w1mag > (w1mag - 16.31) * 1.))
            | (rmag - w1mag > 3.4)
        )  # double sliding cuts and high-z extension
    else:
        raise ValueError

    #################################################################

    mask_lrg &= mask_quality
    if np.sum(mask_lrg)==0:
        return None

    cat = cat[mask_lrg]

    return cat


cat = []

for field in ['north', 'south']:

    sweep_dir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0'.format(field)

    sweep_path_list = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))
    # sweep_fn_list = [os.path.basename(sweep_path_list[ii]) for ii in range(len(sweep_path_list))]

    # start multiple worker processes
    with Pool(processes=n_processes) as pool:
        res = pool.map(get_lrgs, sweep_path_list)

    # res = []
    # for index, sweep_path in enumerate(sweep_path_list):
    #     print(index, sweep_path)
    #     res.append(get_lrgs(sweep_path))

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

cat.write('/global/cfs/cdirs/desi/users/rongpu/test/desi_lrg_targets.fits')

print('All done!', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))

