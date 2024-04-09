# Combine the randoms to speed up the read_pixel_nexp and read_pixel_bitmask steps

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
from multiprocessing import Pool

cat = []
randoms_sizes = []
for randoms_index in range(10):
    tmp = Table(fitsio.read(f'/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-{randoms_index}.fits', columns=['RA', 'DEC', 'BRICKID']))
    randoms_sizes.append(len(tmp))
    cat.append(tmp)
cat = vstack(cat)
print(len(cat))

print(randoms_sizes)
cat.write('/pscratch/sd/r/rongpu/lrg_xcorr/dr10_photoz/randoms/dr9_randoms.fits')

# Run:
# python read_pixel_nexp.py --dr 10 --input /pscratch/sd/r/rongpu/lrg_xcorr/dr10_photoz/randoms/dr9_randoms.fits --output /pscratch/sd/r/rongpu/lrg_xcorr/dr10_photoz/randoms/dr9_randoms_dr10_nexp.fits

randoms_sizes = [51738616, 51738616, 51738616, 51738616, 51738616, 51738616, 51738616, 51738616, 51738616, 51738616]
cumsum = np.cumsum(randoms_sizes)

cat = Table(fitsio.read('/pscratch/sd/r/rongpu/lrg_xcorr/dr10_photoz/randoms/dr9_randoms_dr10_nexp.fits'))
cat.rename_columns(['PIXEL_NOBS_G', 'PIXEL_NOBS_R', 'PIXEL_NOBS_I', 'PIXEL_NOBS_Z'], ['DR10_NOBS_G', 'DR10_NOBS_R', 'DR10_NOBS_I', 'DR10_NOBS_Z'])
assert len(cat)==cumsum[-1]
# for randoms_index in range(10):
#     if randoms_index==0:
#         tmp = cat[0:cumsum[randoms_index]].copy()
#     else:
#         tmp = cat[cumsum[randoms_index-1]:cumsum[randoms_index]].copy()
#     tmp.write(f'/global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr10_photoz/randoms/randoms-1-{randoms_index}-dr10_nexp.fits.gz')


def write_randoms_catalogs(randoms_index):
    if randoms_index==0:
        tmp = cat[0:cumsum[randoms_index]].copy()
    else:
        tmp = cat[cumsum[randoms_index-1]:cumsum[randoms_index]].copy()
    tmp.write(f'/global/cfs/cdirs/desicollab/users/rongpu/data/lrg_xcorr/catalogs/dr10_photoz/randoms/randoms-1-{randoms_index}-dr10_nexp.fits.gz')


n_process = 10
with Pool(processes=n_process) as pool:
    res = pool.map(write_randoms_catalogs, np.arange(10))
