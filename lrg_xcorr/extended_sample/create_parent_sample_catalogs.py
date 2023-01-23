from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

from multiprocessing import Pool


basic_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'PHOTSYS', 'TARGETID', 'SV1_DESI_TARGET', 'SV1_BGS_TARGET']
photom_columns = ['MORPHTYPE', 'EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'FLUX_W1', 'FLUX_W2', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG']

target_columns_all = basic_columns + photom_columns
target_columns_all = list(set(target_columns_all))  # unique columns

target_path_list = glob.glob('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/targets/sv1/resolve/dark/sv1targets-dark-hp-*.fits')

target_bit = 6  # LRG_SV_IR


def read_target_files(target_path):
    tmp = fitsio.read(target_path, columns=['SV1_DESI_TARGET'])
    mask = ((tmp["SV1_DESI_TARGET"] & (2**target_bit))!=0)
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
print(len(cat))

# Require zfibermag<22.0
zfibermag = 22.5 - 2.5*np.log10(cat['FIBERFLUX_Z']) - 1.211 * cat['EBV']
mask = zfibermag<22.0
cat = cat[mask]
print(len(cat))

cat_basic = cat[basic_columns].copy()
cat_photom = cat[photom_columns].copy()

cat_basic = vstack(cat_basic, join_type='exact')
cat_photom = vstack(cat_photom, join_type='exact')
print(len(cat_basic))

cat_basic.write('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/more/dr9_extended_lrg_0.49.0_basic.fits')
cat_photom.write('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/more/dr9_extended_lrg_0.49.0_photom.fits')


