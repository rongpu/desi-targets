from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp


min_nobs = 1
maskbits = [1, 8, 9, 11, 12, 13]
custom_mask_name = ''
mask_str = ''.join([str(tmp) for tmp in maskbits])
if custom_mask_name!='':
    mask_str += '_' + custom_mask_name

fn = '/pscratch/sd/r/rongpu/ebv/g_limited_sample.fits'
cat = Table(fitsio.read(fn))

nside = 128
npix = hp.nside2npix(nside)

################################################################################

pixmap_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/resolve/combined/old'

maps_dict = {}
for field in ['north', 'south']:
    maps = Table.read(os.path.join(pixmap_dir, 'pixmap_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, mask_str)))
    maps = maps[maps['FRACAREA']>0]
    # Load stellar density map
    stardens = np.load('/global/cfs/cdirs/desi/users/rongpu/useful/healpix_maps/pixweight-dr7.1-0.22.0_stardens_{}_ring.npy'.format(nside))
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

################################################################################

mask = cat['gmag_sel'].copy()
pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'][mask], cat['DEC'][mask], lonlat=True)
pix_unique, pix_count = np.unique(pix_allobj, return_counts=True)
hp_table = Table()
hp_table['HPXPIXEL'] = pix_unique
# hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_unique, nest=False, lonlat=True)
hp_table['n_obj'] = pix_count
hp_table = join(hp_table, maps, join_type='inner', keys='HPXPIXEL')
hp_table.write('/pscratch/sd/r/rongpu/ebv/density_map_gmag_selection_{}.fits'.format(nside))

mask = cat['gfiber_sel'].copy()
pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'][mask], cat['DEC'][mask], lonlat=True)
pix_unique, pix_count = np.unique(pix_allobj, return_counts=True)
hp_table = Table()
hp_table['HPXPIXEL'] = pix_unique
# hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_unique, nest=False, lonlat=True)
hp_table['n_obj'] = pix_count
hp_table = join(hp_table, maps, join_type='inner', keys='HPXPIXEL')
hp_table.write('/pscratch/sd/r/rongpu/ebv/density_map_gfiber_selection_{}.fits'.format(nside))

