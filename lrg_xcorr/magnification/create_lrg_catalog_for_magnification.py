# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python create_lrg_catalog_for_magnification.py north

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio
from multiprocessing import Pool

from desitarget.targets import encode_targetid

n_processes = 64

field = str(sys.argv[1])
print(field)

sweep_dir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0'.format(field)
pz_dir = '/global/cfs/cdirs/desi/users/rongpu/ls_dr9.0_desi_photoz/pz/'+field

sweep_columns = ['RELEASE', 'BRICKID', 'OBJID', 'TYPE', 'RA', 'DEC', 'DCHISQ', 'EBV',
'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2',
'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2',
'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2',
'FIBERFLUX_Z', 'FIBERTOTFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG',
'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS',
'SHAPE_R', 'SHAPE_R_IVAR', 'SHAPE_E1', 'SHAPE_E2', 'SERSIC']


def get_lrgs(sweep_fn):

    sweep_path = os.path.join(sweep_dir, sweep_fn)
    pz_path = os.path.join(pz_dir, sweep_fn[:-5]+'-pz.fits')

    cat = Table(fitsio.read(sweep_path, columns=sweep_columns))

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

    mask_lrg_all = np.full(len(cat), False)

    for magnification in [0.99, 1., 1.01]:
        gmag = 22.5 - 2.5 * np.log10((cat['FLUX_G'] * magnification / cat['MW_TRANSMISSION_G']).clip(1e-7))
        # ADM safe as these fluxes are set to > 0 in notinLRG_mask.
        rmag = 22.5 - 2.5 * np.log10((cat['FLUX_R'] * magnification / cat['MW_TRANSMISSION_R']).clip(1e-7))
        zmag = 22.5 - 2.5 * np.log10((cat['FLUX_Z'] * magnification / cat['MW_TRANSMISSION_Z']).clip(1e-7))
        w1mag = 22.5 - 2.5 * np.log10((cat['FLUX_W1'] * magnification / cat['MW_TRANSMISSION_W1']).clip(1e-7))

        # To enable adjustments to the size effect, no change in zfibermag if it's getting fainter
        if magnification<=1:
            zfibermag = 22.5 - 2.5 * np.log10((cat['FIBERFLUX_Z'] / cat['MW_TRANSMISSION_Z']).clip(1e-7))
        else:
            zfibermag = 22.5 - 2.5 * np.log10((cat['FIBERFLUX_Z'] * magnification / cat['MW_TRANSMISSION_Z']).clip(1e-7))

        mask_lrg = np.full(len(cat), True)

        if field=='south':
            mask_lrg &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
            mask_lrg &= zfibermag < 21.6                   # faint limit
            mask_lrg &= (gmag - w1mag > 2.9) | (rmag - w1mag > 1.8)  # low-z cuts
            mask_lrg &= (
                ((rmag - w1mag > (w1mag - 17.14) * 1.8)
                 & (rmag - w1mag > (w1mag - 16.33) * 1.))
                | (rmag - w1mag > 3.3)
            )  # double sliding cuts and high-z extension
        else:
            mask_lrg &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
            mask_lrg &= zfibermag < 21.61                   # faint limit
            mask_lrg &= (gmag - w1mag > 2.97) | (rmag - w1mag > 1.8)  # low-z cuts
            mask_lrg &= (
                ((rmag - w1mag > (w1mag - 17.13) * 1.83)
                 & (rmag - w1mag > (w1mag - 16.31) * 1.))
                | (rmag - w1mag > 3.4)
            )  # double sliding cuts and high-z extension

        mask_lrg_all |= mask_lrg
        print('magnification = {}; '.format(magnification), np.sum(mask_lrg_all))

    mask_lrg_all &= mask_quality
    if np.sum(mask_lrg_all)==0:
        return None

    cat = cat[mask_lrg_all]

    idx = np.where(mask_lrg_all)[0]
    pz = Table(fitsio.read(pz_path, rows=idx))

    pz.remove_columns(['OBJID', 'BRICKID', 'RELEASE'])
    cat = hstack([cat, pz], join_type='exact')

    return cat


print('Start!')

time_start = time.time()

sweep_all_path = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))
sweep_fn_list = [os.path.basename(sweep_all_path[ii]) for ii in range(len(sweep_all_path))]

# start multiple worker processes
with Pool(processes=n_processes) as pool:
    res = pool.map(get_lrgs, sweep_fn_list)

# res = []
# for index, sweep_fn in enumerate(sweep_fn_list):
#     print(index, sweep_fn)
#     res.append(get_lrgs(sweep_fn))

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

cat = vstack(res)

cat['TARGETID'] = encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE'])

cat.write('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/lrg_magnification_{}.fits'.format(field), overwrite=True)

print('All done!', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))

