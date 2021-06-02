# Example:
# python compute_density_variations-sv3.py LRG south

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

# from multiprocessing import Pool
import healpy as hp


target_class, field = str(sys.argv[1]), str(sys.argv[2])
target_class = target_class.upper()
field = field.lower()

min_nobs = 1

target_bits = {'LRG': 0, 'ELG': 1, 'QSO': 2, 'BGS_ANY': 60, 'BGS_BRIGHT': 1}
maskbits_dict = {'LRG': [1, 8, 9, 11, 12, 13], 'ELG': [1, 11, 12, 13], 'QSO': [1, 8, 9, 11, 12, 13], 'BGS_ANY': [1, 13], 'BGS_BRIGHT': [1, 13]}

nsides = [64, 128, 256, 512]

target_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS']

if 'BGS' in target_class:
    target_dir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.57.0/targets/sv3/resolve/bright'
else:
    target_dir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.57.0/targets/sv3/resolve/dark'

output_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/density_maps/0.57.0/resolve'

target_bit = target_bits[target_class]
maskbits = maskbits_dict[target_class]

if field=='south':
    photsys = 'S'
else:
    photsys = 'N'


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
    print(target_class, field)

    time_start = time.time()

    # Load targets
    target_path_list = glob.glob(os.path.join(target_dir, 'sv3targets-*.fits'))
    cat = []
    for target_path in target_path_list:
        # print(target_path)
        if target_class!='BGS_BRIGHT':
            tmp = fitsio.read(target_path, columns=['SV3_DESI_TARGET', 'PHOTSYS'])
            mask = ((tmp["SV3_DESI_TARGET"] & (2**target_bit))!=0) & (tmp['PHOTSYS']==photsys)
        else:
            tmp = fitsio.read(target_path, columns=['SV3_BGS_TARGET', 'PHOTSYS'])
            mask = ((tmp["SV3_BGS_TARGET"] & (2**target_bit))!=0) & (tmp['PHOTSYS']==photsys)
        idx = np.where(mask)[0]
        if len(idx)==0:
            continue
        # print(len(idx)/len(tmp), len(idx), len(tmp))
        cat.append(Table(fitsio.read(target_path, columns=target_columns, rows=idx)))
    cat = vstack(cat)

    print('Loading complete!')

    mask = apply_mask(cat, min_nobs, maskbits)
    cat = cat[mask]

    for nside in nsides:
        npix = hp.nside2npix(nside)
        pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'], cat['DEC'], lonlat=True)
        pix_unique, pix_count = np.unique(pix_allobj, return_counts=True)
        hp_table = get_systeamtics(pix_unique)
        hp_table['n_targets'] = pix_count
        hp_table.write(os.path.join(output_dir, 'density_map_sv3_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(target_class.lower(), field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))

    print(time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
