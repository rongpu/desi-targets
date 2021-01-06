# Create DR9 SV BGS catalog for Alexie

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

target_dir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.47.0/targets/sv1/resolve/bright'

target_class = 'BGS_ANY'

# field = 'south'
# if field=='south':
#     photsys = 'S'
# else:
#     photsys = 'N'

target_bits = {'LRG': 0, 'ELG': 1, 'QSO': 2, 'BGS_ANY': 60}
target_bit = target_bits[target_class]

target_columns = ['RA', 'DEC', 'PHOTSYS', 'SV1_DESI_TARGET', 'SV1_BGS_TARGET']

# bgs_bits = {'BGS_FAINT': 0, 'BGS_BRIGHT': 1, 'BGS_WISE': 2, 'BGS_FAINT_HIP': 3}
sv_bgs_bits = {'BGS_FAINT': 0, 'BGS_BRIGHT': 1, 'BGS_FAINT_EXT': 2, 'BGS_LOWQ': 3, 'BGS_FIBMAG': 4}

# Load targets
target_path_list = glob.glob(os.path.join(target_dir, 'sv1targets-bright-hp-*.fits'))
cat = []
for target_path in target_path_list:
    print(target_path)
    # tmp = fitsio.read(target_path, columns=['SV1_DESI_TARGET', 'PHOTSYS'])
    # mask = ((tmp["SV1_DESI_TARGET"] & (2**target_bit))!=0) & (tmp['PHOTSYS']==photsys)
    tmp = fitsio.read(target_path, columns=['SV1_DESI_TARGET'])
    mask = ((tmp["SV1_DESI_TARGET"] & (2**target_bit))!=0)
    idx = np.where(mask)[0]
    if len(idx)==0:
        continue
    print(len(idx)/len(tmp), len(idx), len(tmp))
    cat.append(Table(fitsio.read(target_path, columns=target_columns, rows=idx)))
cat = vstack(cat)

for bgs_type in sv_bgs_bits.keys():
    bgs_bit = sv_bgs_bits[bgs_type]
    mask = (cat["SV1_BGS_TARGET"] & (2**bgs_bit))!=0
    print(bgs_type, np.sum(mask))
    cat[bgs_type] = np.array(mask, dtype=np.int16)

columns_to_keep = ['RA', 'DEC', 'BGS_FAINT', 'BGS_BRIGHT', 'BGS_FAINT_EXT', 'BGS_LOWQ', 'BGS_FIBMAG']
cat = cat[columns_to_keep]

cat.write('/global/cscratch1/sd/rongpu/share/dr9_lrg/dr9_sv_bgs_radec_only.fits')
