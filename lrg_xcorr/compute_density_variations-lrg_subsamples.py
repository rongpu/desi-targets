from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

import healpy as hp

min_nobs = 2

nsides = [64, 128, 256]

output_dir = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/density_maps/main_lrg'


def get_systematics(pix_list):

    hp_table = Table()
    hp_table['HPXPIXEL'] = pix_list
    hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_list, nest=False, lonlat=True)

    # if target_class=='LRG' or target_class=='QSO':
    #     hp_columns = ['NOBS_W1', 'NOBS_W2']
    #     arr = np.zeros([len(pix_list), len(hp_columns)])
    #     hp_table = hstack([hp_table, Table(arr, names=hp_columns)])
    #     for index in range(len(pix_list)):
    #         mask = pix_allobj==pix_list[index]
    #         hp_table['NOBS_W1'][index] = np.mean(cat['NOBS_W1'][mask])
    #         hp_table['NOBS_W2'][index] = np.mean(cat['NOBS_W2'][mask])

    return hp_table


if __name__ == '__main__':

    print('Start!')

    time_start = time.time()

    # Load LRG catalog
    cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/dr9_lrg_pzbins_20230509.fits'))

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

            for nside in nsides:
                output_path = os.path.join(output_dir, 'density_map_lrg_pz_bin_{}_{}_nside_{}_minobs_{}.fits'.format(bin_index, field, nside, min_nobs))
                if os.path.isfile(output_path):
                    continue
                npix = hp.nside2npix(nside)
                pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'][mask], cat['DEC'][mask], lonlat=True)
                pix_unique, pix_count = np.unique(pix_allobj, return_counts=True)
                hp_table = get_systematics(pix_unique)
                hp_table['n_targets'] = pix_count
                hp_table.write(output_path)

    print('Done!', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
