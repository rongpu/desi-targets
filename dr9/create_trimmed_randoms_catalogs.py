# No mask applied
# Example:
# python create_trimmed_randoms_catalogs.py south

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack
import fitsio

min_nobs = 1
randoms_path = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-0.fits'

output_dir = '/global/cscratch1/sd/rongpu/target/catalogs/dr9.0/0.49.0'

randoms_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z',
                   'GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z', 'PSFDEPTH_W1', 'PSFDEPTH_W2', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z',
                   'APFLUX_G', 'APFLUX_R', 'APFLUX_Z', 'APFLUX_IVAR_G', 'APFLUX_IVAR_R', 'APFLUX_IVAR_Z', 'MASKBITS',
                   'WISEMASK_W1', 'WISEMASK_W2', 'EBV', 'PHOTSYS']

# field: "north" or "south"
field = str(sys.argv[1])
print(field)

if field=='south':
    photsys = 'S'
elif field=='north':
    photsys = 'N'

output_path = os.path.join(output_dir, 'dr9_randoms_{}_0.49.0.fits'.format(field))
if os.path.isfile(output_path):
    sys.exit('File already exist: '+output_path)

randoms = fitsio.read(randoms_path, columns=randoms_columns)

mask = (randoms['PHOTSYS']==photsys)
randoms = randoms[mask]

mask = (randoms['NOBS_G']>=min_nobs) & (randoms['NOBS_R']>=min_nobs) & (randoms['NOBS_Z']>=min_nobs)
print('NOBS cut: {} ({:.2f}%) removed'.format(np.sum(~mask), np.sum(~mask)/len(mask)*100))
randoms = randoms[mask]
print(len(randoms))

header_dict = {'DENSITY':2500}
fitsio.write(output_path, randoms, clobber=False, header=header_dict)
