# Compute linear regression weight maps using healpix maps

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp

sys.path.append(os.path.expanduser('~/git/desi-targets/dr9/density_variations/'))
from compute_weights import get_weights


def write_output(maps):

    print(len(maps))

    area = np.sum(maps['FRACAREA'])*pix_area
    print('Area = {:.1f} sq deg'.format(area))

    mask = maps['FRACAREA']>min_pix_frac
    maps = maps[mask]

    maps['density'] = maps['n_targets'] / (pix_area * maps['FRACAREA'])
    maps['density_predict'] = maps['n_targets_predict'] / (pix_area * maps['FRACAREA'])

    if combine_north_and_south:
        output_path = os.path.join(output_dir, '{}/linear_weights_{}_sv3_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(weights_ver, weights_ver, target_class.lower(), 'combined', nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])+lrgmask_str))
    else:
        output_path = os.path.join(output_dir, '{}/linear_weights_{}_sv3_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(weights_ver, weights_ver, target_class.lower(), field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])+lrgmask_str))
    print(output_path)

    if os.path.isfile(output_path):
        return None

    if not os.path.isdir(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))

    maps.write(output_path)


combine_north_and_south = True

min_nobs = 1
# maskbits_dict = {'LRG': [1, 8, 9, 11, 12, 13], 'ELG': [1, 11, 12, 13], 'QSO': [1, 8, 9, 11, 12, 13], 'BGS_ANY': [1, 13], 'BGS_BRIGHT': [1, 13]}
maskbits_dict = {'LRG': [], 'ELG': [1, 11, 12, 13], 'QSO': [1, 8, 9, 11, 12, 13], 'BGS_ANY': [1, 13], 'BGS_BRIGHT': [1, 13]}
target_ver_str = '0.57.0'
resolve = 'resolve'

randoms_counts_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/{}/counts'.format(resolve)
randoms_systematics_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/{}/systematics'.format(resolve)
target_densities_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/density_maps/{}/{}'.format(target_ver_str, resolve)

output_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/density_maps/{}/{}/linear_weights'.format(target_ver_str, resolve)

min_pix_frac = 0.  # minimum fraction of pixel area to be used
                   # set to 0. to remove pixels that have targets but no randoms

# for index in range(len(xnames)):
#     print(xnames[index], xlabels[index])

weights_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/linear_weights'
# weights_dir = '/Users/rongpu/Documents/Data/desi_targets/dr9.0/imaging_sys/linear_weights'
weights_fn_dict = {'LRG': 'sv3_v0.2/sv3_lrg_linear_coeffs_v0.2.yaml',
                   'BGS_ANY': 'sv3_v0.1/sv3_bgs_any_linear_coeffs_v0.1.yaml',
                   'BGS_BRIGHT': 'sv3_v0.1/sv3_bgs_bright_linear_coeffs_v0.1.yaml',
                   'ELG': 'sv3_v0.1/sv3_elg_linear_coeffs_v0.1.yaml',
                   'QSO': 'sv3_v0.1/sv3_qso_linear_coeffs_separate_des_v0.1.yaml',
                   }

for target_class in ['BGS_ANY', 'BGS_BRIGHT', 'ELG', 'QSO', 'LRG']:

    if target_class=='LRG':
        lrgmask_str = '_lrgmask_v1'
    else:
        lrgmask_str = ''

    print(target_class)
    target_class = target_class.lower()

    maskbits = maskbits_dict[target_class.upper()]

    weights_fn = weights_fn_dict[target_class.upper()]
    weights_path = os.path.join(weights_dir, weights_fn)
    if '_separate_des_' in os.path.basename(weights_fn):
        separate_des = True
    else:
        separate_des = False
    print("Separate DES weights:", separate_des)

    if target_class.upper()=='LRG':
        weights_ver = 'v0.2'
    else:
        weights_ver = 'v0.1'

    for nside in [64, 128, 256, 512]:

        npix = hp.nside2npix(nside)
        pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
        print(nside, 'Healpix size = {:.5f} sq deg'.format(pix_area))

        field = 'combined'

        for field in ['north', 'south']:

            if field=='north' or field=='BASS+MzLS':
                photsys = 'N'
            else:
                photsys = 'S'

            density = Table.read(os.path.join(target_densities_dir, 'density_map_sv3_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(target_class, field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])+lrgmask_str)))
            maps = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])+lrgmask_str)))
            maps = maps[maps['n_randoms']>0]
            maps1 = Table.read(os.path.join(randoms_systematics_dir, 'systematics_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])+lrgmask_str)))
            maps1.remove_columns(['RA', 'DEC'])
            maps = join(maps, maps1, join_type='inner', keys='HPXPIXEL')
            maps = join(maps, density[['HPXPIXEL', 'n_targets']], join_type='outer', keys='HPXPIXEL').filled(0)

            # # Load stellar density map
            # stardens = np.load('/Users/rongpu/Documents/Data/desi_lrg_selection/dr7/healpix_maps/pixweight-dr7.1-0.22.0_stardens_{}_ring.npy'.format(nside))
            # maps['stardens'] = stardens[maps['HPXPIXEL']]
            # maps['stardens_log'] = np.log10(maps['stardens'])

            maps['PHOTSYS'] = photsys

            density_predict = get_weights(maps, weights_path, separate_des=separate_des)
            maps['n_targets_predict'] = density_predict * (pix_area * maps['FRACAREA'])

            if not combine_north_and_south:
                write_output(maps)
            else:
                if field=='north':
                    maps_north = maps.copy()
                else:
                    maps_south = maps.copy()

        if combine_north_and_south:

            ########## Combine the two maps; proper handling of overlapping pixels ##########

            pix_overlap = np.intersect1d(maps_north['HPXPIXEL'], maps_south['HPXPIXEL'])
            mask = np.in1d(maps_north['HPXPIXEL'], pix_overlap)
            maps_overlap_north = maps_north[mask]
            maps_north = maps_north[~mask]
            mask = np.in1d(maps_south['HPXPIXEL'], pix_overlap)
            maps_overlap_south = maps_south[mask]
            maps_south = maps_south[~mask]

            maps_overlap_north.sort('HPXPIXEL')
            maps_overlap_south.sort('HPXPIXEL')

            maps_overlap = maps_overlap_south.copy()
            maps_overlap['n_targets'] = maps_overlap_north['n_targets'] + maps_overlap_south['n_targets']
            maps_overlap['n_targets_predict'] = maps_overlap_north['n_targets_predict'] + maps_overlap_south['n_targets_predict']
            maps_overlap['FRACAREA'] = maps_overlap_north['FRACAREA'] + maps_overlap_south['FRACAREA']

            maps = vstack([maps_north, maps_south, maps_overlap])

            ######################################################################

            write_output(maps)

