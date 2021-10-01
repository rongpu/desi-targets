# version 1.1

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

# from multiprocessing import Pool
import healpy as hp

min_nobs = 2

nsides = [64, 128, 256, 512]

output_dir = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/density_maps/1.0.0/resolve/v1.1'

# maskbits = [1, 8, 9, 11, 12, 13]

maskbits = []
apply_lrgmask = True
if apply_lrgmask:
    lrgmask_str = '_lrgmask_v1'
else:
    lrgmask_str = ''


def apply_mask(cat, min_nobs, maskbits):

    mask = (cat['NOBS_G']>=min_nobs) & (cat['NOBS_R']>=min_nobs) & (cat['NOBS_Z']>=min_nobs)

    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    # print(np.sum(~mask_clean)/len(mask_clean))

    mask &= mask_clean

    return mask


def get_systeamtics(pix_list):

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

    for field in ['north', 'south']:

        if field=='south':
            photsys = 'S'
        elif field=='north':
            photsys = 'N'

        # Load LRG catalog
        min_nobs_cat = 1
        # cat_path = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/main_lrg_minobs_{}_maskbits_{}_20210723.fits'.format(min_nobs_cat, ''.join([str(tmp) for tmp in maskbits]))
        cat_path = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/main_lrg_minobs_{}_20210913.fits'.format(min_nobs_cat)
        cat = Table(fitsio.read(cat_path))

        mask = cat['PHOTSYS']==photsys
        cat = cat[mask]

        mask = cat['lrg_mask']==0
        cat = cat[mask]

        mask = apply_mask(cat, min_nobs, maskbits)
        cat = cat[mask]
        
        ############################## photo-z bins ##############################

        for bin_index in range(1, 5):  # 4 bins

            mask = cat['pz_bin']==bin_index

            for nside in nsides:
                # output_path = os.path.join(output_dir, 'density_map_lrg_pz_bin_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(bin_index, field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])))
                output_path = os.path.join(output_dir, 'density_map_lrg_pz_bin_{}_{}_nside_{}_minobs_{}.fits'.format(bin_index, field, nside, min_nobs))
                if os.path.isfile(output_path):
                    continue
                npix = hp.nside2npix(nside)
                pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'][mask], cat['DEC'][mask], lonlat=True)
                pix_unique, pix_count = np.unique(pix_allobj, return_counts=True)
                hp_table = get_systeamtics(pix_unique)
                hp_table['n_targets'] = pix_count
                hp_table.write(output_path)

    print('Done!', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
