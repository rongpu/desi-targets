# Subsample version 1.0

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

from multiprocessing import Pool
import healpy as hp

import yaml

n_processes = 32

min_nobs = 2

nsides = [64, 128, 256, 512]
# nsides = [64]

weights_path = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/imaging_weights/v1.0/main_lrg_linear_coeffs_pz.yaml'

output_dir = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/density_maps/1.0.0/resolve/v1.0/linear_weights'

maskbits = [1, 8, 9, 11, 12, 13]


def apply_mask(cat, min_nobs, maskbits):

    mask = (cat['NOBS_G']>=min_nobs) & (cat['NOBS_R']>=min_nobs) & (cat['NOBS_Z']>=min_nobs)

    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    # print(np.sum(~mask_clean)/len(mask_clean))

    mask &= mask_clean

    return mask


def get_weighted_counts(pix_idx):

    pix_list = pix_unique[pix_idx]

    hp_table = Table()
    hp_table['HPXPIXEL'] = pix_list
    hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_list, nest=False, lonlat=True)
    hp_table['n_targets'] = 0.

    for index in np.arange(len(pix_idx)):

        idx = pixorder[pixcnts[pix_idx[index]]:pixcnts[pix_idx[index]+1]]
        hp_table['n_targets'][index] = np.sum(weights[idx])

    return hp_table


if __name__ == '__main__':

    print('Start!')

    time_start = time.time()

    for field in ['north', 'south']:

        # Load weights
        with open(weights_path, "r") as f:
            linear_coeffs = yaml.safe_load(f)

        # Load LRG catalog
        min_nobs_cat = 1
        cat_path = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/main_lrg_minobs_{}_maskbits_{}_20210723.fits'.format(min_nobs_cat, ''.join([str(tmp) for tmp in maskbits]))
        cat_more_path = cat_path.replace('.fits', '_more.fits')
        cat = Table(fitsio.read(cat_path))
        cat_more = Table(fitsio.read(cat_more_path))
        cat_more.remove_columns(['TARGETID', 'EBV'])
        cat = hstack([cat, cat_more])

        if field=='south':
            photsys = 'S'
        elif field=='north':
            photsys = 'N'
        mask = cat['PHOTSYS']==photsys
        cat = cat[mask]

        mask = apply_mask(cat, min_nobs, maskbits)
        cat = cat[mask]

        cat_all = cat.copy()
        
        ############################## photo-z bins ##############################

        for bin_index in range(1, 5):  # 4 bins

            mask = cat_all['pz_bin']==bin_index
            cat = cat_all[mask].copy()

            xnames_fit = list(linear_coeffs['south_bin_1'].keys())
            xnames_fit.remove('intercept')
            # Assign zero weights to objects with invalid imaging properties
            # (their fraction should be negligibly small)
            mask_bad = np.full(len(cat), False)
            for col in xnames_fit:
                mask_bad |= ~np.isfinite(cat[col])
            if np.sum(mask_bad)!=0:
                print('{} invalid objects'.format(np.sum(mask_bad)))

            weights = np.zeros(len(cat))
            bin_str = '{}_bin_{}'.format(field, bin_index)

            # create array of coefficients, with the first coefficient being the intercept
            coeffs = np.array([linear_coeffs[bin_str]['intercept']]+[linear_coeffs[bin_str][xname] for xname in xnames_fit])

            data = np.column_stack([cat[~mask_bad][xname] for xname in xnames_fit])
            # create 2-D array of imaging properties, with the first columns being unity
            data1 = np.insert(data, 0, 1., axis=1)
            # wt = coeff0 + coeff1 * rand['EBV'] + coeff2 * rand['PSFSIZE_G'] + ...
            weights[~mask_bad] = 1/np.dot(coeffs, data1.T)  # 1/predicted_density as weights for objects

            for nside in nsides:
                output_path = os.path.join(output_dir, 'density_map_lrg_pz_bin_{}_{}_nside_{}_minobs_{}_maskbits_{}-lw.fits'.format(bin_index, field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])))
                if os.path.isfile(output_path):
                    continue
                npix = hp.nside2npix(nside)

                pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'], cat['DEC'], lonlat=True)
                pix_unique, pixcnts = np.unique(pix_allobj, return_counts=True)

                pixcnts = np.insert(pixcnts, 0, 0)
                pixcnts = np.cumsum(pixcnts)

                pixorder = np.argsort(pix_allobj)

                # split among the Cori processors
                pix_idx_split = np.array_split(np.arange(len(pix_unique)), n_processes)

                # start multiple worker processes
                with Pool(processes=n_processes) as pool:
                    res = pool.map(get_weighted_counts, pix_idx_split)

                hp_table = vstack(res)
                hp_table.sort('HPXPIXEL')

                hp_table.write(output_path)

    print('Done!', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
