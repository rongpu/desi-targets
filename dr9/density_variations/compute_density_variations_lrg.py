from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

# from multiprocessing import Pool
import healpy as hp

field = str(sys.argv[1])
field = field.lower()

min_nobs = 1
maskbits = sorted([1, 8, 9, 11, 12, 13])

n_processes = 32

nsides = [64, 128, 256]

target_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'NOBS_W1']

target_path = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/unofficial/sv3_lrg_{}.fits'.format(field)

output_dir = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/unofficial/density_maps'

# hp_columns = ['NOBS_W1', 'NOBS_W2']
hp_columns = ['NOBS_W1']


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
    hp_table['hp_idx'] = pix_list
    hp_table['ra'], hp_table['dec'] = hp.pixelfunc.pix2ang(nside, pix_list, nest=False, lonlat=True)

    arr = np.zeros([len(pix_list), len(hp_columns)])
    hp_table = hstack([hp_table, Table(arr, names=hp_columns)])

    # Use mean for NOBS; median for other properties
    for index in range(len(pix_list)):

        mask = pix_allobj==pix_list[index]
        hp_table['NOBS_W1'][index] = np.mean(cat['NOBS_W1'][mask])
        # hp_table['NOBS_W2'][index] = np.mean(cat['NOBS_W2'][mask])

    return hp_table


if __name__ == '__main__':

    print('Start!')

    time_start = time.time()

    cat = Table(fitsio.read(target_path))
    print('Loading complete!')

    mask = apply_mask(cat, min_nobs, maskbits)
    cat = cat[mask]

    # 800 targets/sq.deg. selection
    for nside in nsides:
        npix = hp.nside2npix(nside)
        pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'], cat['DEC'], lonlat=True)
        pix_unique, pix_count = np.unique(pix_allobj, return_counts=True)
        hp_table = get_systeamtics(pix_unique)
        hp_table['count'] = pix_count
        hp_table.write(os.path.join(output_dir, 'density_map_sv3_lrg_all_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))

    # 600 targets/sq.deg. selection
    mask = cat['lrg_lowdens'].copy()
    cat = cat[mask]
    for nside in nsides:
        npix = hp.nside2npix(nside)
        pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'], cat['DEC'], lonlat=True)
        pix_unique, pix_count = np.unique(pix_allobj, return_counts=True)
        hp_table = get_systeamtics(pix_unique)
        hp_table['count'] = pix_count
        hp_table.write(os.path.join(output_dir, 'density_map_sv3_lrg_lowdens_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))

    print(time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))

