# Example:
# salloc -N 1 -C cpu -q interactive -t 4:00:00
# python add_pixel_level_values.py LRG

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack
import fitsio

from multiprocessing import Pool

target_pix_dir = '/global/cfs/cdirs/desi/survey/catalogs/extra_target_data/1.1.1'
output_dir = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.1.1/resolve'

# pixel_columns = ['TARGETID', 'DESI_TARGET', 'PIXEL_NOBS_G', 'PIXEL_NOBS_R', 'PIXEL_NOBS_Z', 'PIXEL_PSFDEPTH_G', 'PIXEL_PSFDEPTH_R', 'PIXEL_PSFDEPTH_Z', 'PIXEL_GALDEPTH_G', 'PIXEL_GALDEPTH_R', 'PIXEL_GALDEPTH_Z', 'PIXEL_PSFDEPTH_W1', 'PIXEL_PSFDEPTH_W2', 'PIXEL_PSFSIZE_G', 'PIXEL_PSFSIZE_R', 'PIXEL_PSFSIZE_Z', 'PIXEL_APFLUX_G', 'PIXEL_APFLUX_R', 'PIXEL_APFLUX_Z', 'PIXEL_APFLUX_IVAR_G', 'PIXEL_APFLUX_IVAR_R', 'PIXEL_APFLUX_IVAR_Z', 'PIXEL_MASKBITS', 'PIXEL_WISEMASK_W1', 'PIXEL_WISEMASK_W2']
pixel_columns = ['TARGETID', 'PIXEL_NOBS_G', 'PIXEL_NOBS_R', 'PIXEL_NOBS_Z']

# target_class: "LRG", "ELG", "QSO" or "BGS_ANY"
target_class = str(sys.argv[1])
target_class = target_class.upper()

print(target_class)

if target_class=='BGS_ANY':
    target_path_list = glob.glob(os.path.join(target_pix_dir, 'bright', '*.fits'))
else:
    target_path_list = glob.glob(os.path.join(target_pix_dir, 'dark', '*.fits'))

cat_basic_path = os.path.join(output_dir, 'dr9_{}_1.1.1_basic.fits'.format(target_class.lower()))
cat_pixel_path = os.path.join(output_dir, 'dr9_{}_1.1.1_pixel.fits'.format(target_class.lower()))

if os.path.isfile(cat_pixel_path):
    sys.exit('File already exist: '+cat_pixel_path)

cat_basic = Table(fitsio.read(cat_basic_path, columns=['RA', 'DEC', 'TARGETID', 'PHOTSYS']))


def read_target_files(target_path):
    tmp = fitsio.read(target_path, columns=['TARGETID'])
    mask = np.in1d(tmp['TARGETID'], cat_basic['TARGETID'])
    idx = np.where(mask)[0]
    if len(idx)==0:
        return None
    cat = Table(fitsio.read(target_path, columns=pixel_columns, rows=idx))
    return cat


n_processes = 128
with Pool(processes=n_processes) as pool:
    res = pool.map(read_target_files, target_path_list, chunksize=1)

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

cat = vstack(res, join_type='exact')

# Here matching cat to cat_basic
t1_reverse_sort = np.array(cat_basic['TARGETID']).argsort().argsort()
cat = cat[np.argsort(cat['TARGETID'])[t1_reverse_sort]]
if not np.all(cat['TARGETID']==cat_basic['TARGETID']):
    raise ValueError('different targetid')
cat.remove_column('TARGETID')

cat.write(cat_pixel_path)
