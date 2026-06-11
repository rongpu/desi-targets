# Create SV3 LRG catalog with no masking and "noresolve"

from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
import gc
from multiprocessing import Pool

from desitarget.targets import encode_targetid


field = str(sys.argv[1])
field = field.lower()
n_processes = 32

sweep_dir = '/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0'.format(field)
output_path = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/unofficial/sv3_lrg_{}.fits'.format(field)

sweep_all_path = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))
sweep_fn_list = [os.path.basename(sweep_all_path[ii]) for ii in range(len(sweep_all_path))]

columns = ['OBJID', 'BRICKID', 'RELEASE',
           'TYPE', 'RA', 'DEC', 'EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1',
           'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'NOBS_W1',
           'FRACFLUX_G', 'FRACFLUX_R', 'FRACFLUX_Z', 'FRACFLUX_W1', 'FRACIN_G', 'FRACIN_R', 'FRACIN_Z', 'WISEMASK_W1',
           'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z', 'PSFDEPTH_W1',
           'SHAPE_R', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG',
           'GAIA_ASTROMETRIC_EXCESS_NOISE', 'MASKBITS', 'FITBITS', 'SERSIC']


def trim_sweep(sweep_index):

    sweep_fn = sweep_fn_list[sweep_index]
    print(sweep_fn)

    cat = fitsio.read(os.path.join(sweep_dir, sweep_fn),
                      columns=['FIBERFLUX_Z', 'FLUX_IVAR_Z', 'MW_TRANSMISSION_Z'])
    cat = Table(cat)
    # Apply fibermag cut zfibermag<21.8
    mask = (cat['FIBERFLUX_Z']>0) & (cat['FLUX_IVAR_Z']>0)
    mask &= (cat['FIBERFLUX_Z']/cat['MW_TRANSMISSION_Z'] > 10**(0.4*(22.5-21.8)))
    if np.sum(mask)==0:
        return None
    idx = np.where(mask)[0]

    cat = fitsio.read(os.path.join(sweep_dir, sweep_fn), columns=columns, rows=idx)
    cat = Table(cat)
    # cat['TYPE'] = cat['TYPE'].astype(str)

    lrg_quality = (cat['FLUX_IVAR_R'] > 0) & (cat['FLUX_R'] > 0)   # ADM quality in r.
    lrg_quality &= (cat['FLUX_IVAR_Z'] > 0) & (cat['FLUX_Z'] > 0) & (cat['FIBERFLUX_Z'] > 0)   # ADM quality in z.
    lrg_quality &= (cat['FLUX_IVAR_W1'] > 0) & (cat['FLUX_W1'] > 0)  # ADM quality in W1.

    lrg_quality &= (cat['GAIA_PHOT_G_MEAN_MAG'] == 0) | (cat['GAIA_PHOT_G_MEAN_MAG'] > 18)  # remove bright GAIA sources

    # ADM remove stars with zfibertot < 17.5 that are missing from GAIA (no extinction correction)
    lrg_quality &= cat['FIBERTOTFLUX_Z'] < 10**(-0.4*(17.5-22.5))

    # ADM observed in every band.
    lrg_quality &= (cat['NOBS_G']>0) & (cat['NOBS_R']>0) & (cat['NOBS_Z']>0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gmag = 22.5 - 2.5*np.log10(cat['FLUX_G']/cat['MW_TRANSMISSION_G'])
        rmag = 22.5 - 2.5*np.log10(cat['FLUX_R']/cat['MW_TRANSMISSION_R'])
        zmag = 22.5 - 2.5*np.log10(cat['FLUX_Z']/cat['MW_TRANSMISSION_Z'])
        w1mag = 22.5 - 2.5*np.log10(cat['FLUX_W1']/cat['MW_TRANSMISSION_W1'])
        rfibermag = 22.5 - 2.5*np.log10(cat['FIBERFLUX_R']/cat['MW_TRANSMISSION_R'])
        zfibermag = 22.5 - 2.5*np.log10(cat['FIBERFLUX_Z']/cat['MW_TRANSMISSION_Z'])

        gmag[~np.isfinite(gmag)] = 100.
        rmag[~np.isfinite(rmag)] = 100.
        zmag[~np.isfinite(zmag)] = 100.
        w1mag[~np.isfinite(w1mag)] = 100.
        rfibermag[~np.isfinite(rfibermag)] = 100.
        zfibermag[~np.isfinite(zfibermag)] = 100.

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # South
        if field=='south':

            # 800 targets/sq.deg. selection

            lrg_mask = lrg_quality.copy()
            lrg_mask &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
            lrg_mask &= (zfibermag < 21.7)                   # faint limit

            lrg_mask_sliding = rmag - w1mag > (w1mag - 17.26) * 1.8  # sliding IR cut
            lrg_mask_sliding &= rmag - w1mag > (w1mag - 16.36) * 1.  # low-z sliding IR cut
            lrg_mask_sliding |= rmag - w1mag > 3.29
            lrg_mask &= lrg_mask_sliding

            mask_lowz = (gmag - rmag > 1.3) & ((gmag - rmag) > -1.55 * (rmag - w1mag) + 3.13)
            mask_lowz |= (rmag - w1mag > 1.8)
            lrg_mask &= mask_lowz

            # 600 targets/sq.deg. selection

            lrg_mask_lowdens = lrg_quality.copy()
            lrg_mask_lowdens &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
            lrg_mask_lowdens &= (zfibermag < 21.7)                   # faint limit

            lrg_mask_sliding = rmag - w1mag > (w1mag - 17.07) * 1.8  # sliding IR cut
            lrg_mask_sliding &= rmag - w1mag > (w1mag - 16.17) * 1.  # low-z sliding IR cut
            lrg_mask_sliding |= rmag - w1mag > 3.39
            lrg_mask_lowdens &= lrg_mask_sliding

            mask_lowz = (gmag - rmag > 1.3) & ((gmag - rmag) > -1.55 * (rmag - w1mag) + 3.13)
            mask_lowz |= (rmag - w1mag > 1.8)
            lrg_mask_lowdens &= mask_lowz

        # North
        elif field=='north':

            # 800 targets/sq.deg. selection

            lrg_mask = lrg_quality.copy()
            lrg_mask &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
            lrg_mask &= (zfibermag < 21.72)                   # faint limit

            lrg_mask_sliding = rmag - w1mag > (w1mag - 17.24) * 1.83  # sliding IR cut
            lrg_mask_sliding &= rmag - w1mag > (w1mag - 16.33) * 1.   # low-z sliding IR cut
            lrg_mask_sliding |= rmag - w1mag > 3.39
            lrg_mask &= lrg_mask_sliding

            mask_lowz = (gmag - rmag > 1.34) & ((gmag - rmag) > -1.55 * (rmag - w1mag) + 3.23)
            mask_lowz |= (rmag - w1mag > 1.8)
            lrg_mask &= mask_lowz

            # 600 targets/sq.deg. selection

            lrg_mask_lowdens = lrg_quality.copy()
            lrg_mask_lowdens &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
            lrg_mask_lowdens &= (zfibermag < 21.72)                   # faint limit

            lrg_mask_sliding = rmag - w1mag > (w1mag - 17.05) * 1.83  # sliding IR cut
            lrg_mask_sliding &= rmag - w1mag > (w1mag - 16.14) * 1.   # low-z sliding IR cut
            lrg_mask_sliding |= rmag - w1mag > 3.49
            lrg_mask_lowdens &= lrg_mask_sliding

            mask_lowz = (gmag - rmag > 1.34) & ((gmag - rmag) > -1.55 * (rmag - w1mag) + 3.23)
            mask_lowz |= (rmag - w1mag > 1.8)
            lrg_mask_lowdens &= mask_lowz

        else:
            raise ValueError('0')

    if np.sum(lrg_mask)>0:
        cat = cat[lrg_mask]
        if np.sum(lrg_mask_lowdens & lrg_mask)!=np.sum(lrg_mask_lowdens):
            raise ValueError('1')
        lrg_mask_lowdens = lrg_mask_lowdens[lrg_mask]
        cat['lrg_lowdens'] = np.full(len(cat), False)
        cat['lrg_lowdens'][lrg_mask_lowdens] = True
    else:
        return None

    cat['TARGETID'] = encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE'])

    # clear cache
    gc.collect()

    return cat


if __name__ == '__main__':

    print('Start!')

    time_start = time.time()

    # start multiple worker processes
    with Pool(processes=n_processes) as pool:
        res = pool.map(trim_sweep, range(len(sweep_fn_list)))

    # Remove None elements from the list
    for index in range(len(res)-1, -1, -1):
        if res[index] is None:
            res.pop(index)

    cat_stack = vstack(res)

    # # "Fix" bug in vstack that creates excessively long strings
    # cat_stack['TYPE'] = cat_stack['TYPE'].astype('a4')

    print('Final combined catalog:', len(cat_stack))

    cat_stack.write(output_path)

    print(time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
