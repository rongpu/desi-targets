from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp

sys.path.append(os.path.expanduser('~/git/desi-examples/misc/misc'))
from healpix_util import downsize_hp_map


nside_in = 512
nside_out = 128
bad_sky_mask = False

min_nobs = 1
maskbits = []
custom_mask_name = 'elgmask_v1'
mask_str = ''.join([str(tmp) for tmp in maskbits])
if custom_mask_name!='':
    mask_str += '_' + custom_mask_name

fn = '/pscratch/sd/r/rongpu/ebv/alternative_elg_targets/alternative_elg_targets.fits'
elgmask_fn = '/pscratch/sd/r/rongpu/ebv/alternative_elg_targets/alternative_elg_targets_elgmask_v1.fits.gz'
cat = Table(fitsio.read(fn))
elgmask = Table(fitsio.read(elgmask_fn))
cat = hstack([cat, elgmask], join_type='exact')
mask = cat['elg_mask']==0
cat = cat[mask]

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

pixmap_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/resolve/combined'

maps_dict = {}
for field in ['north', 'south']:
    maps = Table.read(os.path.join(pixmap_dir, 'pixmap_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside_in, min_nobs, mask_str)))
    maps = maps[maps['FRACAREA']>0]
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
maps_overlap['FRACAREA'] = maps_overlap_north['FRACAREA'] + maps_overlap_south['FRACAREA']

maps = vstack([maps_north, maps_south, maps_overlap])

######################################################################

if bad_sky_mask:
    mask = ~np.in1d(maps['HPXPIXEL'], bad_pixels)
    print('Bad sky', np.sum(mask)/len(mask))
    maps = maps[mask]

print('maps', len(maps))

######################################################################

stats_dict = {'FRACAREA': np.sum, 'n_targets': np.sum, 'EBV': np.average, 'GALDEPTH_G': np.average, 'GALDEPTH_R': np.average, 'GALDEPTH_Z': np.average, 'PSFDEPTH_G': np.average, 'PSFDEPTH_R': np.average, 'PSFDEPTH_Z': np.average, 'PSFDEPTH_W1': np.average, 'PSFDEPTH_W2': np.average, 'PSFSIZE_G': np.average, 'PSFSIZE_R': np.average, 'PSFSIZE_Z': np.average, 'STARDENS': np.average}
columns = ['HPXPIXEL', 'FRACAREA', 'n_targets', 'EBV', 'GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z', 'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'PSFDEPTH_W1', 'PSFDEPTH_W2', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 'STARDENS']
n_processes = 128

mask = cat['elg_original'].copy()
output_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/v0/targets/maps/density_map_alternative_elgs_original_{}.fits'.format(nside_out)
if bad_sky_mask:
    output_fn = output_fn.replace('.fits', '_badskymask.fits')
pix_allobj = hp.pixelfunc.ang2pix(nside_in, cat['RA'][mask], cat['DEC'][mask], lonlat=True)
pix_unique, pix_count = np.unique(pix_allobj, return_counts=True)
hp_table = Table()
hp_table['HPXPIXEL'] = pix_unique
# hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_unique, nest=False, lonlat=True)
hp_table['n_targets'] = pix_count
hp_table = join(maps, hp_table, keys='HPXPIXEL', join_type='left').filled(0)
hp_table = downsize_hp_map(nside_in, nside_out, hp_table[columns], stats_dict=stats_dict, weights=hp_table['FRACAREA'], n_processes=n_processes)
hp_table.rename_column('FRACAREA', 'FRACAREA_IN')
pix_area_in = hp.pixelfunc.nside2pixarea(nside_in, degrees=True)
pix_area_out = hp.pixelfunc.nside2pixarea(nside_out, degrees=True)
hp_table['FRACAREA'] = hp_table['FRACAREA_IN']*pix_area_in/pix_area_out
hp_table.write(output_fn, overwrite=True)

mask = cat['elg_gmag'].copy()
output_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/v0/targets/maps/density_map_alternative_elgs_gmag_{}.fits'.format(nside_out)
if bad_sky_mask:
    output_fn = output_fn.replace('.fits', '_badskymask.fits')
pix_allobj = hp.pixelfunc.ang2pix(nside_in, cat['RA'][mask], cat['DEC'][mask], lonlat=True)
pix_unique, pix_count = np.unique(pix_allobj, return_counts=True)
hp_table = Table()
hp_table['HPXPIXEL'] = pix_unique
# hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_unique, nest=False, lonlat=True)
hp_table['n_targets'] = pix_count
hp_table = join(maps, hp_table, keys='HPXPIXEL', join_type='left').filled(0)
hp_table = downsize_hp_map(nside_in, nside_out, hp_table[columns], stats_dict=stats_dict, weights=hp_table['FRACAREA'], n_processes=n_processes)
hp_table.rename_column('FRACAREA', 'FRACAREA_IN')
pix_area_in = hp.pixelfunc.nside2pixarea(nside_in, degrees=True)
pix_area_out = hp.pixelfunc.nside2pixarea(nside_out, degrees=True)
hp_table['FRACAREA'] = hp_table['FRACAREA_IN']*pix_area_in/pix_area_out
hp_table.write(output_fn, overwrite=True)

mask = cat['elg_brighter'].copy()
output_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/v0/targets/maps/density_map_alternative_elgs_brighter_{}.fits'.format(nside_out)
if bad_sky_mask:
    output_fn = output_fn.replace('.fits', '_badskymask.fits')
pix_allobj = hp.pixelfunc.ang2pix(nside_in, cat['RA'][mask], cat['DEC'][mask], lonlat=True)
pix_unique, pix_count = np.unique(pix_allobj, return_counts=True)
hp_table = Table()
hp_table['HPXPIXEL'] = pix_unique
# hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_unique, nest=False, lonlat=True)
hp_table['n_targets'] = pix_count
hp_table = join(maps, hp_table, keys='HPXPIXEL', join_type='left').filled(0)
hp_table = downsize_hp_map(nside_in, nside_out, hp_table[columns], stats_dict=stats_dict, weights=hp_table['FRACAREA'], n_processes=n_processes)
hp_table.rename_column('FRACAREA', 'FRACAREA_IN')
pix_area_in = hp.pixelfunc.nside2pixarea(nside_in, degrees=True)
pix_area_out = hp.pixelfunc.nside2pixarea(nside_out, degrees=True)
hp_table['FRACAREA'] = hp_table['FRACAREA_IN']*pix_area_in/pix_area_out
hp_table.write(output_fn, overwrite=True)

mask = cat['elg_gmag_brighter'].copy()
output_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/v0/targets/maps/density_map_alternative_elgs_gmag_brighter_{}.fits'.format(nside_out)
if bad_sky_mask:
    output_fn = output_fn.replace('.fits', '_badskymask.fits')
pix_allobj = hp.pixelfunc.ang2pix(nside_in, cat['RA'][mask], cat['DEC'][mask], lonlat=True)
pix_unique, pix_count = np.unique(pix_allobj, return_counts=True)
hp_table = Table()
hp_table['HPXPIXEL'] = pix_unique
# hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_unique, nest=False, lonlat=True)
hp_table['n_targets'] = pix_count
hp_table = join(maps, hp_table, keys='HPXPIXEL', join_type='left').filled(0)
hp_table = downsize_hp_map(nside_in, nside_out, hp_table[columns], stats_dict=stats_dict, weights=hp_table['FRACAREA'], n_processes=n_processes)
hp_table.rename_column('FRACAREA', 'FRACAREA_IN')
pix_area_in = hp.pixelfunc.nside2pixarea(nside_in, degrees=True)
pix_area_out = hp.pixelfunc.nside2pixarea(nside_out, degrees=True)
hp_table['FRACAREA'] = hp_table['FRACAREA_IN']*pix_area_in/pix_area_out
hp_table.write(output_fn, overwrite=True)
