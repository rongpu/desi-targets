from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

unwise_bit = 0
ts_bit = 4

input_dir = '/global/cscratch1/sd/rongpu/desi/lrg_pixel_bitmask/dev'
output_dir = '/global/cscratch1/sd/rongpu/desi/lrg_pixel_bitmask/v1'

field = 'south'

unwise_maskbits = [0, 1, 2, 3, 4, 6, 7]  # all except the HALO bit
ts_maskbits = [1, 12, 13]  # DESI targeting mask bits


def get_combined_bitmask(brick_index):

    brickname = str(bricks['brickname'][brick_index])

    output_path = os.path.join(output_dir, '{}/coadd/{}/{}/{}-lrgmask.fits.gz'.format(field, brickname[:3], brickname, brickname))
    if os.path.isfile(output_path):
        return None

    print(output_path)

    img_fn = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/coadd/{}/{}/legacysurvey-{}-maskbits.fits.fz'.format(field, brickname[:3], brickname, brickname)

    lrgmask = np.full((3600, 3600), 0, dtype=np.int16)

    # unWISE maskbits
    wisem1 = fitsio.read(img_fn, ext='WISEM1')
    mask_unwise = np.full((3600, 3600), False)
    for bit in unwise_maskbits:
        mask_unwise |= (wisem1 & 2**bit)>0
    lrgmask[mask_unwise] += 2**unwise_bit

    # targeting maskbits
    maskbits = fitsio.read(img_fn, ext='MASKBITS')
    mask_ts = np.full(maskbits.shape, False)
    for bit in ts_maskbits:
        mask_ts |= (maskbits & 2**bit)>0
    lrgmask[mask_ts] += 2**ts_bit

    gaiam_path = os.path.join(input_dir, '{}/coadd/{}/{}/{}-wisemask.fits.gz'.format(field, brickname[:3], brickname, brickname))
    wisem_path = os.path.join(input_dir, '{}/coadd/{}/{}/{}-gaiamask.fits.gz'.format(field, brickname[:3], brickname, brickname))

    gaiam = fitsio.read(gaiam_path)
    wisem = fitsio.read(wisem_path)

    lrgmask += gaiam
    lrgmask += wisem

    if not os.path.isdir(os.path.dirname(output_path)):
        try:
            os.makedirs(os.path.dirname(output_path))
        except:
            pass
    fitsio.write(output_path, lrgmask, compress='GZIP')

    return lrgmask


bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/survey-bricks-dr9-{}.fits.gz'.format(field, field)))
print(len(bricks))

mask = bricks['brickname']=='1092p320'
brick_index = np.where(mask)[0][0]

lrgmask = get_combined_bitmask(brick_index)

