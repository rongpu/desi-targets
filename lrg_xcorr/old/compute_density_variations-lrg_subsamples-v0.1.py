# Subsample version 0.1

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

# from multiprocessing import Pool
import healpy as hp

min_nobs = 1

nsides = [64, 128, 256, 512]

basic_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS']
photom_columns = ['EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z',
                  'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'FLUX_W1',
                  'MW_TRANSMISSION_W1', 'FIBERFLUX_Z']

output_dir = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/density_maps/1.0.0/resolve'

maskbits = [1, 8, 9, 11, 12, 13]


def apply_mask(cat, min_nobs, maskbits):

    mask = (cat['NOBS_G']>=min_nobs) & (cat['NOBS_R']>=min_nobs) & (cat['NOBS_Z']>=min_nobs)

    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    # print(np.sum(~mask_clean)/len(mask_clean))

    mask &= mask_clean

    return mask


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

    for field in ['north', 'south']:

        # Load targets
        cat0 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve/dr9_lrg_{}_1.0.0_basic.fits'.format(field), columns=basic_columns))
        cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve/dr9_lrg_{}_1.0.0_photom.fits'.format(field), columns=photom_columns))
        pz = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve/dr9_lrg_{}_1.0.0_pz.fits'.format(field)))
        cat = hstack([cat0, cat, pz], join_type='exact')
        print(len(cat))

        print('Loading complete!')

        mask = apply_mask(cat, min_nobs, maskbits)
        cat = cat[mask]

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gmag = 22.5 - 2.5*np.log10(cat['FLUX_G']/cat['MW_TRANSMISSION_G'])
            rmag = 22.5 - 2.5*np.log10(cat['FLUX_R']/cat['MW_TRANSMISSION_R'])
            zmag = 22.5 - 2.5*np.log10(cat['FLUX_Z']/cat['MW_TRANSMISSION_Z'])
            w1mag = 22.5 - 2.5*np.log10(cat['FLUX_W1']/cat['MW_TRANSMISSION_W1'])
            zfibermag = 22.5 - 2.5*np.log10(cat['FIBERFLUX_Z']/cat['MW_TRANSMISSION_Z'])

        ############################## r-W1 bins ##############################

        if field=='south':
            rw1_cuts = [2.24, 2.82, 3.3, 3.54]
        else:
            rw1_cuts = [2.32, 2.92, 3.4, 3.65]

        for bin_index in range(-1, len(rw1_cuts)):

            if bin_index==-1:
                rw1_min, rw1_max = -np.inf, rw1_cuts[0]
            elif bin_index==len(rw1_cuts)-1:
                rw1_min, rw1_max = rw1_cuts[bin_index], np.inf
            else:
                rw1_min, rw1_max = rw1_cuts[bin_index], rw1_cuts[bin_index+1]
            mask = (rmag-w1mag>rw1_min) & (rmag-w1mag<rw1_max)

            for nside in nsides:
                output_path = os.path.join(output_dir, 'density_map_lrg_rw1_bin_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(bin_index+2, field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])))
                if os.path.isfile(output_path):
                    continue
                npix = hp.nside2npix(nside)
                pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'][mask], cat['DEC'][mask], lonlat=True)
                pix_unique, pix_count = np.unique(pix_allobj, return_counts=True)
                hp_table = get_systematics(pix_unique)
                hp_table['n_targets'] = pix_count
                hp_table.write(output_path)

        ############################## photo-z bins ##############################

        if field=='south':
            pz_cuts = [0.540, 0.683, 0.810, 0.890]
        else:
            pz_cuts = [0.552, 0.691, 0.812, 0.885]

        for bin_index in range(-1, len(pz_cuts)):

            if bin_index==-1:
                pz_min, pz_max = 0, pz_cuts[0]
            elif bin_index==len(pz_cuts)-1:
                pz_min, pz_max = pz_cuts[bin_index], np.inf
            else:
                pz_min, pz_max = pz_cuts[bin_index], pz_cuts[bin_index+1]
            mask = (cat['Z_PHOT_MEDIAN']>pz_min) & (cat['Z_PHOT_MEDIAN']<pz_max)

            for nside in nsides:
                output_path = os.path.join(output_dir, 'density_map_lrg_pz_bin_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(bin_index+2, field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])))
                if os.path.isfile(output_path):
                    continue
                npix = hp.nside2npix(nside)
                pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'][mask], cat['DEC'][mask], lonlat=True)
                pix_unique, pix_count = np.unique(pix_allobj, return_counts=True)
                hp_table = get_systematics(pix_unique)
                hp_table['n_targets'] = pix_count
                hp_table.write(output_path)

    print('Done!', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
