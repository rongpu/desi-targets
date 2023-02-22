from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp
from multiprocessing import Pool

sys.path.append(os.path.expanduser('~/git/desi-examples/misc/misc'))
from healpix_util import downsize_hp_map


target_class = 'ELG_LOP'
nside_in = 1024
nside_out = 128
bad_sky_mask = True

######################################################################

if bad_sky_mask:
    bad_pixels = Table.read('/global/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/bad_pixels_v1_512_ring.fits')['HPXPIXEL']
    if nside_in!=512:
        npix = hp.nside2npix(512)
        tmp = np.zeros(npix)
        tmp[bad_pixels] = 1.
        tmp = hp.ud_grade(tmp, nside_in, order_in='RING', order_out='RING')
        bad_pixels = np.where(tmp!=0)[0]

######################################################################


target_ver_str = '1.1.1'
pixmap_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/resolve/combined'
target_densities_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/density_maps/{}/resolve'.format(target_ver_str)

min_nobs = 1

maskbits_dict = {'LRG': [], 'ELG': [], 'ELG_LOP': [], 'QSO': [1, 8, 9, 11, 12, 13], 'BGS_ANY': [1, 13], 'BGS_BRIGHT': [1, 13]}
custom_mask_dict = {'LRG': 'lrgmask_v1.1', 'ELG': 'elgmask_v1', 'ELG_LOP': 'elgmask_v1', 'QSO': '', 'BGS_ANY': '', 'BGS_BRIGHT': ''}

maskbits = maskbits_dict[target_class]
custom_mask_name = custom_mask_dict[target_class]

mask_str = ''.join([str(tmp) for tmp in maskbits])
if custom_mask_name!='':
    mask_str += '_' + custom_mask_name

target_class = target_class.lower()
maps_dict = {}

for field in ['north', 'south']:

    density = Table.read(os.path.join(target_densities_dir, 'density_map_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(target_class, field, nside_in, min_nobs, mask_str)))
    maps = Table.read(os.path.join(pixmap_dir, 'pixmap_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside_in, min_nobs, mask_str)))
    maps = maps[maps['FRACAREA']>0]
    maps = join(maps, density[['HPXPIXEL', 'n_targets']], join_type='inner', keys='HPXPIXEL').filled(0)

    # Load stellar density map
    stardens = np.load('/global/cfs/cdirs/desi/users/rongpu/useful/healpix_maps/pixweight-dr7.1-0.22.0_stardens_{}_ring.npy'.format(nside_in))
    maps['stardens'] = stardens[maps['HPXPIXEL']]
    maps['stardens_log'] = np.log10(maps['stardens'])

    maps_dict[field] = maps.copy()

maps_north = maps_dict['north']
maps_south = maps_dict['south']

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
maps_overlap['FRACAREA'] = maps_overlap_north['FRACAREA'] + maps_overlap_south['FRACAREA']

maps = vstack([maps_north, maps_south, maps_overlap])

######################################################################

fn_hi = '/global/cfs/cdirs/desi/users/rongpu/useful/lenz_hi/ebv_lhd_equatorial.hpx.fits'
ebv = Table(fitsio.read(fn_hi))
ebv['HPXPIXEL'] = np.arange(hp.nside2npix(nside_in))
ebv.rename_column('EBV', 'EBV_HI')
mask = np.isfinite(ebv['EBV_HI'])
ebv = ebv[mask]
mask = np.in1d(maps['HPXPIXEL'], ebv['HPXPIXEL'])
print('HI EBV map', np.sum(mask)/len(mask))
maps = maps[mask]
maps = join(maps, ebv, keys='HPXPIXEL', join_type='inner')

######################################################################

if bad_sky_mask:
    mask = ~np.in1d(maps['HPXPIXEL'], bad_pixels)
    print('Bad sky', np.sum(mask)/len(mask))
    maps = maps[mask]

print('maps', len(maps))

######################################################################

new_targets = Table(fitsio.read('/pscratch/sd/r/rongpu/ebv/count_map_{}_{}.fits'.format(target_class.lower(), nside_in)))
new_targets = new_targets[['HPXPIXEL', 'n_obj']]
new_targets.rename_column('n_obj', 'n_targets_new')

maps = join(maps, new_targets, keys='HPXPIXEL', join_type='left').filled(0)

stats_dict = {'FRACAREA': np.sum, 'n_targets': np.sum, 'n_targets_new': np.sum, 'EBV': np.average, 'EBV_HI': np.average, 'GALDEPTH_G': np.average, 'GALDEPTH_R': np.average, 'GALDEPTH_Z': np.average, 'PSFDEPTH_G': np.average, 'PSFDEPTH_R': np.average, 'PSFDEPTH_Z': np.average, 'PSFDEPTH_W1': np.average, 'PSFDEPTH_W2': np.average, 'PSFSIZE_G': np.average, 'PSFSIZE_R': np.average, 'PSFSIZE_Z': np.average, 'STARDENS': np.average}
columns = ['HPXPIXEL', 'FRACAREA', 'n_targets', 'n_targets_new', 'EBV', 'EBV_HI', 'GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z', 'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'PSFDEPTH_W1', 'PSFDEPTH_W2', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 'STARDENS']

n_processes = 128
maps_ds = downsize_hp_map(nside_in, nside_out, maps[columns], stats_dict=stats_dict, weights=maps['FRACAREA'], n_processes=n_processes)
maps_ds.rename_column('FRACAREA', 'FRACAREA_IN')

pix_area_in = hp.pixelfunc.nside2pixarea(nside_in, degrees=True)
pix_area_out = hp.pixelfunc.nside2pixarea(nside_out, degrees=True)
maps_ds['FRACAREA'] = maps_ds['FRACAREA_IN']*pix_area_in/pix_area_out

maps_ds['density'] = maps_ds['n_targets'] / (pix_area_out * maps_ds['FRACAREA'])
maps_ds['density_new'] = maps_ds['n_targets_new'] / (pix_area_out * maps_ds['FRACAREA'])

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    maps_ds['galdepth_gmag'] = -2.5*(np.log10((5/np.sqrt(maps_ds['GALDEPTH_G'])))-9)
    maps_ds['galdepth_rmag'] = -2.5*(np.log10((5/np.sqrt(maps_ds['GALDEPTH_R'])))-9)
    maps_ds['galdepth_zmag'] = -2.5*(np.log10((5/np.sqrt(maps_ds['GALDEPTH_Z'])))-9)
    maps_ds['psfdepth_gmag'] = -2.5*(np.log10((5/np.sqrt(maps_ds['PSFDEPTH_G'])))-9)
    maps_ds['psfdepth_rmag'] = -2.5*(np.log10((5/np.sqrt(maps_ds['PSFDEPTH_R'])))-9)
    maps_ds['psfdepth_zmag'] = -2.5*(np.log10((5/np.sqrt(maps_ds['PSFDEPTH_Z'])))-9)
    maps_ds['psfdepth_w1mag'] = -2.5*(np.log10((5/np.sqrt(maps_ds['PSFDEPTH_W1'])))-9)
    maps_ds['psfdepth_w2mag'] = -2.5*(np.log10((5/np.sqrt(maps_ds['PSFDEPTH_W2'])))-9)
    # maps_ds['galdepth_gmag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps_ds['GALDEPTH_G'])))-9) - 3.214*maps_ds['EBV']
    # maps_ds['galdepth_rmag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps_ds['GALDEPTH_R'])))-9) - 2.165*maps_ds['EBV']
    # maps_ds['galdepth_zmag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps_ds['GALDEPTH_Z'])))-9) - 1.211*maps_ds['EBV']
    # maps_ds['psfdepth_gmag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps_ds['PSFDEPTH_G'])))-9) - 3.214*maps_ds['EBV']
    # maps_ds['psfdepth_rmag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps_ds['PSFDEPTH_R'])))-9) - 2.165*maps_ds['EBV']
    # maps_ds['psfdepth_zmag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps_ds['PSFDEPTH_Z'])))-9) - 1.211*maps_ds['EBV']
    # maps_ds['psfdepth_w1mag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps_ds['PSFDEPTH_W1'])))-9) - 0.184*maps_ds['EBV']
    # maps_ds['psfdepth_w2mag_ebv'] = -2.5*(np.log10((5/np.sqrt(maps_ds['PSFDEPTH_W2'])))-9) - 0.113*maps_ds['EBV']

output_fn = '/pscratch/sd/r/rongpu/ebv/density_map_{}_{}.fits'.format(target_class.lower(), nside_out)
if bad_sky_mask:
    output_fn = output_fn.replace('.fits', '_badskymask.fits')
maps_ds.write(output_fn)

