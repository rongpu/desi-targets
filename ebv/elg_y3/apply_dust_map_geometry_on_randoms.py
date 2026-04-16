import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
import healpy as hp

# cat_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/misc/desi_targets_with_desi_ebv_y3/elg_targets-desi_egr-256.fits'
cat_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/misc/desi_targets_with_desi_ebv_y3/elg_targets-desi_egr_erz-256.fits'

filename = os.path.basename(cat_fn)
nside = int(filename.split('-')[-1].split('.')[0])
print(nside)

if filename.split('-')[-2]=='desi_egr':

    dust_map = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars_y3/v0.2/final_maps/desi_dust_gr_{}.fits'.format(nside)))
    hpix = np.array(dust_map['HPXPIXEL'])

elif filename.split('-')[-2]=='desi_egr_erz':

    dust_map_gr = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars_y3/v0.2/final_maps/desi_dust_gr_{}.fits'.format(nside)))
    dust_map_rz = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars_y3/v0.2/final_maps/desi_dust_rz_{}.fits'.format(nside)))
    hpix = np.intersect1d(dust_map_gr['HPXPIXEL'], dust_map_rz['HPXPIXEL'])

else:
    raise ValueError('Something is wrong')


# Load randoms catalog
columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'EBV']
randoms = Table(fitsio.read('/dvs_ro/cfs/cdirs/desi/public/ets/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-0.fits', columns=columns))

# Remove randoms outside the DESI dust map
randoms['HPXPIXEL'] = hp.ang2pix(nside, randoms['RA'], randoms['DEC'], nest=False, lonlat=True)
mask_rand = np.in1d(randoms['HPXPIXEL'], hpix)
print(np.sum(mask_rand)/len(mask_rand))
randoms = randoms[mask_rand]
