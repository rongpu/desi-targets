# Example:
# python create_per_tracer_catalogs.py LRG south

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack
import fitsio

min_nobs = 1
maskbits_dark = [1, 12, 13]  # Default maskbits for LRG/ELG/QSO ["BRIGHT", "GALAXY", "CLUSTER"]
maskbits_bgs = [1, 13]  # Default maskbits for bgs ["BRIGHT", "CLUSTER"]

target_dir = '/global/cfs/cdirs/desi/target/catalogs/dr9/1.0.0/targets/main/resolve'
output_dir = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.0.0/resolve'

basic_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'PHOTSYS', 'TARGETID', 'DESI_TARGET', 'BGS_TARGET']
photom_columns = ['MORPHTYPE', 'EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z',
                  'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'FLUX_W1', 'FLUX_W2',
                  'FLUX_IVAR_W1', 'FLUX_IVAR_W2', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2',
                  'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'SHAPE_R', 'SERSIC', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG']

# target_class: "LRG", "ELG", "QSO" or "BGS_ANY"
# field: "north" or "south"
target_class, field = str(sys.argv[1]), str(sys.argv[2])
target_class = target_class.upper()
field = field.lower()

# The following target bits are the same in both main and SV3
target_bits = {'LRG': 0, 'ELG': 1, 'QSO': 2, 'BGS_ANY': 60}
target_bit = target_bits[target_class]

if field=='south':
    photsys = 'S'
elif field=='north':
    photsys = 'N'

# #####################
# target_class = 'LRG'
# #####################

print(target_class, field)

if target_class=='BGS_ANY':
    target_path_list = glob.glob(os.path.join(target_dir, 'bright', '*.fits'))
    maskbits = maskbits_bgs
else:
    target_path_list = glob.glob(os.path.join(target_dir, 'dark', '*.fits'))
    maskbits = maskbits_dark

cat_basic_path = os.path.join(output_dir, 'dr9_{}_{}_1.0.0_basic.fits'.format(target_class.lower(), field))
cat_photom_path = os.path.join(output_dir, 'dr9_{}_{}_1.0.0_photom.fits'.format(target_class.lower(), field))

if os.path.isfile(cat_basic_path):
    sys.exit('File already exist: '+cat_basic_path)

# ####################################################################################################################################
# target_path_list = target_path_list[::50]
# ####################################################################################################################################

cat_basic = []
cat_photom = []

for index, target_path in enumerate(target_path_list):
    tmp = fitsio.read(target_path, columns=['DESI_TARGET', 'PHOTSYS'])
    mask = ((tmp["DESI_TARGET"] & (2**target_bit))!=0) & (tmp['PHOTSYS']==photsys)
    idx = np.where(mask)[0]
    if len(idx)==0:
        continue
    # print(index, '/', len(target_path_list), len(idx)/len(tmp), len(idx), len(tmp))
    cat_basic.append(Table(fitsio.read(target_path, columns=basic_columns, rows=idx)))
    cat_photom.append(Table(fitsio.read(target_path, columns=photom_columns, rows=idx)))

cat_basic = vstack(cat_basic, join_type='exact')
cat_photom = vstack(cat_photom, join_type='exact')

mask = (cat_basic['NOBS_G']>=min_nobs) & (cat_basic['NOBS_R']>=min_nobs) & (cat_basic['NOBS_Z']>=min_nobs)
print('NOBS cut: {} ({:.2f}%) removed'.format(np.sum(~mask), np.sum(~mask)/len(mask)*100))
cat_basic = cat_basic[mask]
cat_photom = cat_photom[mask]

# Sanity check
mask_clean = np.ones(len(cat_basic), dtype=bool)
for bit in maskbits:
    mask_clean &= (cat_basic['MASKBITS'] & 2**bit)==0
# print('MASKBITS cut: {} ({:.2f}%) removed'.format(np.sum(~mask_clean), np.sum(~mask_clean)/len(mask_clean)*100))
if np.sum(~mask_clean)>0:
    raise ValueError

cat_basic = cat_basic[mask_clean]
cat_photom = cat_photom[mask_clean]
print(len(cat_basic))
print()

cat_basic.write(cat_basic_path)
cat_photom.write(cat_photom_path)
