# Example:
# salloc -N 1 -C cpu -q interactive -t 4:00:00
# python create_lge_catalog.py LGE

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack
import fitsio

from multiprocessing import Pool

target_dir = '/global/cfs/cdirs/desi/target/catalogs/dr9/3.0.0/targets/main/resolve'
output_dir = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/3.0.0/resolve'

basic_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'PHOTSYS', 'TARGETID', 'DESI_TARGET', 'BGS_TARGET']

photom_columns = ['MORPHTYPE', 'EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z',
                  'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'FLUX_W1', 'FLUX_W2',
                  'FLUX_IVAR_W1', 'FLUX_IVAR_W2', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2',
                  'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG']

target_columns_all = basic_columns + photom_columns
target_columns_all = list(set(target_columns_all))  # unique columns

# target_class: "LGE"
target_class = str(sys.argv[1])
target_class = target_class.upper()

# The following target bits are the same in both main and SV3
target_bits = {'LGE': 3}
target_bit = target_bits[target_class]

print(target_class)

if target_class=='BGS_ANY':
    target_path_list = glob.glob(os.path.join(target_dir, 'bright', '*.fits'))
else:
    target_path_list = glob.glob(os.path.join(target_dir, 'dark1b', '*.fits'))

cat_basic_path = os.path.join(output_dir, 'dr9_{}_3.0.0_basic.fits'.format(target_class.lower()))
cat_photom_path = os.path.join(output_dir, 'dr9_{}_3.0.0_photom.fits'.format(target_class.lower()))

if os.path.isfile(cat_basic_path):
    sys.exit('File already exist: '+cat_basic_path)


def read_target_files(target_path):
    tmp = fitsio.read(target_path, columns=['DESI_TARGET'])
    mask = ((tmp["DESI_TARGET"] & (2**target_bit))!=0)
    idx = np.where(mask)[0]
    if len(idx)==0:
        return None
    cat = Table(fitsio.read(target_path, columns=target_columns_all, rows=idx))
    return cat


n_processes = 128
with Pool(processes=n_processes) as pool:
    res = pool.map(read_target_files, target_path_list, chunksize=1)

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

cat = vstack(res, join_type='exact')

cat_basic = cat[basic_columns].copy()
cat_photom = cat[photom_columns].copy()

cat_basic = vstack(cat_basic, join_type='exact')
cat_photom = vstack(cat_photom, join_type='exact')
print(len(cat_basic))

cat_basic.write(cat_basic_path)
cat_photom.write(cat_photom_path)
