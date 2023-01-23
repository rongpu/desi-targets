# srun -N 1 -C cpu -c 128 -t 04:00:00 -L cfs -q interactive python compute_density_variations-extended_subsamples-weighted.py

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

from multiprocessing import Pool
import healpy as hp

import yaml

n_processes = 128

min_nobs = 2

nsides = [64, 128, 256, 512]
# nsides = [64]

include_ebv = False
output_dir = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/density_maps/extended_lrg/linear_weights'

if include_ebv:
    weights_path = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/more/dr9_extended_lrg_0.49.0_pzbins_20230120-weights.fits'
else:
    weights_path = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/more/dr9_extended_lrg_0.49.0_pzbins_20230120-weights_no_ebv.fits'


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

    # Load LRG catalog
    cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/dr9_extended_lrg_0.49.0_pzbins_20230120.fits'))
    cat_weights = Table(fitsio.read(weights_path))
    cat = hstack([cat, cat_weights], join_type='exact')

    mask = (cat['PIXEL_NOBS_G']>=min_nobs) & (cat['PIXEL_NOBS_R']>=min_nobs) & (cat['PIXEL_NOBS_Z']>=min_nobs)
    mask &= cat['lrg_mask']==0
    cat = cat[mask]

    for field in ['north', 'south']:

        if field=='south':
            photsys = 'S'
        elif field=='north':
            photsys = 'N'

        ############################## photo-z bins ##############################

        for bin_index in range(1, 5):  # 4 bins

            mask = cat['PHOTSYS']==photsys
            mask &= cat['pz_bin']==bin_index
            weights = cat['weight'][mask].copy()

            for nside in nsides:
                if include_ebv:
                    output_path = os.path.join(output_dir, 'density_map_extended_lrg_pz_bin_{}_{}_nside_{}_minobs_{}-lw.fits'.format(bin_index, field, nside, min_nobs))
                else:
                    output_path = os.path.join(output_dir, 'density_map_extended_lrg_pz_bin_{}_{}_nside_{}_minobs_{}-lw_no_ebv.fits'.format(bin_index, field, nside, min_nobs))
                if os.path.isfile(output_path):
                    continue
                npix = hp.nside2npix(nside)

                pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'][mask], cat['DEC'][mask], lonlat=True)
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
