# Use the pixel-level NOBS instead of the tractor NOBS
# Example:
# python compute_systematics_maps_from_targets.py LRG south

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

# from multiprocessing import Pool
import healpy as hp
from multiprocessing import Pool


n_processes = 128

target_class, field = str(sys.argv[1]), str(sys.argv[2])
target_class = target_class.upper()
field = field.lower()

if target_class=='BGS_BRIGHT':
    target_class = 'BGS_ANY'
    sub_class = 'BGS_BRIGHT'
elif target_class=='ELG_LOP':
    target_class = 'ELG'
    sub_class = 'ELG_LOP'
else:
    sub_class = target_class

min_nobs = 1

# maskbits_dict = {'LRG': [1, 12, 13], 'ELG': [1, 12, 13], 'ELG_LOP': [1, 12, 13], 'QSO': [1, 12, 13], 'BGS_ANY': [1, 13], 'BGS_BRIGHT': [1, 13]}  # The maskbits in desitarget
# maskbits_dict = {'LRG': [1, 8, 9, 11, 12, 13], 'ELG': [1, 11, 12, 13], 'QSO': [1, 8, 9, 11, 12, 13], 'BGS_ANY': [1, 13], 'BGS_BRIGHT': [1, 13]}
# custom_mask_dict = {'LRG': '', 'ELG': '', 'ELG_LOP': '', 'QSO': '', 'BGS_ANY': '', 'BGS_BRIGHT': ''}

maskbits_dict = {'LRG': [], 'ELG': [], 'ELG_LOP': [], 'QSO': [1, 8, 9, 11, 12, 13], 'BGS_ANY': [1, 13], 'BGS_BRIGHT': [1, 13]}
custom_mask_dict = {'LRG': 'lrgmask_v1.1', 'ELG': 'elgmask_v1', 'ELG_LOP': 'elgmask_v1', 'QSO': '', 'BGS_ANY': '', 'BGS_BRIGHT': ''}

# nsides = [64, 128, 256, 512, 1024]
nsides = [64, 128, 256, 512]

# target_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS']
# pix_columns = ['PIXEL_NOBS_G', 'PIXEL_NOBS_R', 'PIXEL_NOBS_Z']

if 'BGS' in target_class:
    target_dir = '/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/bright'
else:
    target_dir = '/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/dark'

cat_dir = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.1.1/resolve'

output_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/density_maps/1.1.1/resolve'

maskbits = maskbits_dict[target_class]
custom_mask_name = custom_mask_dict[target_class]

mask_str = ''.join([str(tmp) for tmp in maskbits])
if custom_mask_name!='':
    mask_str += '_' + custom_mask_name


if field=='south':
    photsys = 'S'
else:
    photsys = 'N'

hp_columns = ['EBV', 'GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z', 'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'PSFDEPTH_W1', 'PSFDEPTH_W2', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 'NOBS_G', 'NOBS_R', 'NOBS_Z']


def apply_mask(cat, min_nobs, maskbits, custom_mask_name):

    mask = (cat['PIXEL_NOBS_G']>=min_nobs) & (cat['PIXEL_NOBS_R']>=min_nobs) & (cat['PIXEL_NOBS_Z']>=min_nobs)

    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    # print(np.sum(~mask_clean)/len(mask_clean))

    if custom_mask_name!='':
        mask_col = custom_mask_name[: custom_mask_name.find("mask")]+'_mask'
        mask_clean &= cat[mask_col]==0

    mask &= mask_clean

    return mask


def get_systematics(pix_idx):

    pix_list = pix_unique[pix_idx]

    hp_table = Table()
    hp_table['HPXPIXEL'] = pix_list
    hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_list, nest=False, lonlat=True)

    arr = np.zeros([len(pix_list), len(hp_columns)])
    hp_table = hstack([hp_table, Table(arr, names=hp_columns)])

    for index in np.arange(len(pix_idx)):

        idx = pixorder[pixcnts[pix_idx[index]]:pixcnts[pix_idx[index]+1]]

        for hp_column in hp_columns:
            if 'NOBS_' in hp_column:
                hp_table[hp_column][index] = np.mean(cat[hp_column][idx])
            else:
                hp_table[hp_column][index] = np.median(cat[hp_column][idx])

    return hp_table


if __name__ == '__main__':

    print('Start!')
    print(sub_class, field)

    time_start = time.time()

    cat_path = os.path.join(cat_dir, 'dr9_{}_1.1.1_basic.fits'.format(target_class.lower()))
    pix_path = os.path.join(cat_dir, 'dr9_{}_1.1.1_pixel.fits'.format(target_class.lower()))
    photom_path = os.path.join(cat_dir, 'dr9_{}_1.1.1_photom.fits'.format(target_class.lower()))
    sweep2_path = os.path.join(cat_dir, 'dr9_{}_1.1.1_sweep_2.fits'.format(target_class.lower()))

    cat = Table(fitsio.read(cat_path))
    cat_pix = Table(fitsio.read(pix_path))
    photom = Table(fitsio.read(photom_path))
    sweep2 = Table(fitsio.read(sweep2_path))
    cat = hstack([cat, cat_pix, photom, sweep2], join_type='exact')

    photsys_mask = cat['PHOTSYS']==photsys
    cat = cat[photsys_mask]

    if custom_mask_name!='':
        mask_path = os.path.join(cat_dir, 'dr9_{}_1.1.1_{}.fits.gz'.format(target_class.lower(), custom_mask_name))
        cat_mask = Table(fitsio.read(mask_path))[photsys_mask]
        cat = hstack([cat, cat_mask], join_type='exact')

    if sub_class=='BGS_BRIGHT':
        mask = cat['BGS_TARGET'] & 2**1 > 0
        cat = cat[mask]
    elif sub_class=='ELG_LOP':
        mask = cat['DESI_TARGET'] & 2**5 > 0
        cat = cat[mask]

    print('Loading complete!')

    mask = apply_mask(cat, min_nobs, maskbits, custom_mask_name)
    cat = cat[mask]

    for nside in nsides:

        output_path = os.path.join(output_dir, 'systematics_map_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(sub_class.lower(), field, nside, min_nobs, mask_str))
        if os.path.isfile(output_path):
            continue

        npix = hp.nside2npix(nside)

        pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'], cat['DEC'], lonlat=True)
        pix_unique, pixcnts = np.unique(pix_allobj, return_counts=True)
        pix_count = pixcnts.copy()

        pixcnts = np.insert(pixcnts, 0, 0)
        pixcnts = np.cumsum(pixcnts)

        pixorder = np.argsort(pix_allobj)

        # split among the Cori processors
        pix_idx_split = np.array_split(np.arange(len(pix_unique)), n_processes)

        # start multiple worker processes
        with Pool(processes=n_processes) as pool:
            res = pool.map(get_systematics, pix_idx_split)

        hp_table = vstack(res)
        hp_table.sort('HPXPIXEL')

        hp_table['n_targets'] = pix_count
        hp_table.write(output_path)

    print(time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
