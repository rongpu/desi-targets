# Example script for Julien

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

########################################### Load data ###########################################

randoms_counts_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/resolve/counts'
randoms_systematics_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/resolve/systematics'
target_densities_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/density_maps/1.0.0/resolve'

target_class = 'lrg'
maskbits = []
lrgmask_str = '_lrgmask_v1'
min_nobs = 1
min_pix_frac = 0.2  # minimum fraction of pixel area to be used


nside = 64

npix = hp.nside2npix(nside)
pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
print(nside, 'Healpix size = {:.5f} sq deg'.format(pix_area))

for field in ['north', 'south']:

    density = Table.read(os.path.join(target_densities_dir, 'density_map_{}_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(target_class, field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])+lrgmask_str)))
    maps = Table.read(os.path.join(randoms_counts_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])+lrgmask_str)))
    maps = maps[maps['n_randoms']>0]
    maps1 = Table.read(os.path.join(randoms_systematics_dir, 'systematics_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])+lrgmask_str)))
    maps1.remove_columns(['RA', 'DEC'])
    maps = join(maps, maps1, join_type='inner', keys='HPXPIXEL')
    maps = join(maps, density[['HPXPIXEL', 'n_targets']], join_type='outer', keys='HPXPIXEL').filled(0)

    # # Load stellar density map
    # stardens = np.load('/global/cfs/cdirs/desi/users/rongpu/useful/healpix_maps/pixweight-dr7.1-0.22.0_stardens_{}_ring.npy'.format(nside))
    # maps['stardens'] = stardens[maps['HPXPIXEL']]
    # maps['stardens_log'] = np.log10(maps['stardens'])

    if field=='north':
        maps_north = maps.copy()
    else:
        maps_south = maps.copy()

#### Combine the two maps; proper handling of overlapping pixels ###

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
print(len(maps))

##################

area = np.sum(maps['FRACAREA'])*pix_area
print('Area = {:.1f} sq deg'.format(area))

mask = maps['FRACAREA']>min_pix_frac
maps = maps[mask]

maps['density'] = maps['n_targets'] / (pix_area * maps['FRACAREA'])

########################################### Plot density map ###########################################

map_values = np.zeros(npix)
hp_mask = np.zeros(npix, dtype=bool)
map_values[maps['HPXPIXEL']] = maps['density']
hp_mask[maps['HPXPIXEL']] = True
mplot = hp.ma(map_values)
mplot.mask = ~hp_mask

plt.figure(figsize=(12, 8))
hp.mollview(mplot, min=300, max=900, rot=(120, 0, 0), fig=1, xsize=1000, cmap='jet')
plt.savefig('lrg_density_map.png', dpi=200)
plt.show()

