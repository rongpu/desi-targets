# To run on single interactive node:
# python combine_lrg_masks_rerun.py "south 1 0"

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from multiprocessing import Pool
import argparse


unwise_bit = 0
ts_bit = 4
unwise_maskbits = [0, 1, 2, 3, 4, 6, 7]  # all except the HALO bit
ts_maskbits = [1, 12, 13]  # DESI targeting mask bits

input_dir = '/global/cscratch1/sd/rongpu/desi/lrg_pixel_bitmask/dev'
output_dir = '/global/cscratch1/sd/rongpu/desi/lrg_pixel_bitmask/v1'

# field = 'south'

# parser = argparse.ArgumentParser()
# parser.add_argument('field')
# parser.add_argument('n_task')
# parser.add_argument('task_id')
# args = parser.parse_args()
# field = args.field
# n_task = int(args.n_task)
# task_id = int(args.task_id)

parser = argparse.ArgumentParser()
parser.add_argument('args')
args = parser.parse_args()
field, n_task, task_id = args.args.split()
n_task, task_id = int(n_task), int(task_id)

n_processes = 32


################################
debug = False
################################


def get_combined_bitmask(brick_index):

    brickname = str(bricks['brickname'][brick_index])

    output_path = os.path.join(output_dir, '{}/coadd/{}/{}/{}-lrgmask.fits.gz'.format(field, brickname[:3], brickname, brickname))
    if os.path.isfile(output_path):
        return None

    # print(output_path)

    img_fn = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/coadd/{}/{}/legacysurvey-{}-maskbits.fits.fz'.format(field, brickname[:3], brickname, brickname)

    header = fitsio.read_header(img_fn, ext=1)
    header_keywords = ['CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']
    hdict = {}
    for keyword in header_keywords:
        hdict[keyword] = header[keyword]
    if (hdict['CTYPE1']!='RA---TAN') or (hdict['CTYPE2']!='DEC--TAN'):
        raise ValueError

    bitmask = np.full((3600, 3600), 0, dtype=np.uint8)

    # unWISE maskbits
    wisem1 = fitsio.read(img_fn, ext='WISEM1')
    mask_unwise = np.full((3600, 3600), False)
    for bit in unwise_maskbits:
        mask_unwise |= (wisem1 & 2**bit)>0
    bitmask[mask_unwise] += 2**unwise_bit

    # targeting maskbits
    maskbits = fitsio.read(img_fn, ext='MASKBITS')
    mask_ts = np.full(maskbits.shape, False)
    for bit in ts_maskbits:
        mask_ts |= (maskbits & 2**bit)>0
    bitmask[mask_ts] += 2**ts_bit

    gaiam_path = os.path.join(input_dir, '{}/coadd/{}/{}/{}-gaiamask.fits.gz'.format(field, brickname[:3], brickname, brickname))
    wisem_path = os.path.join(input_dir, '{}/coadd/{}/{}/{}-wisemask.fits.gz'.format(field, brickname[:3], brickname, brickname))

    gaiam = fitsio.read(gaiam_path).astype(np.uint8)
    wisem = fitsio.read(wisem_path).astype(np.uint8)

    bitmask += gaiam
    bitmask += wisem

    if not os.path.isdir(os.path.dirname(output_path)):
        try:
            os.makedirs(os.path.dirname(output_path))
        except:
            pass
    fitsio.write(output_path, bitmask, compress='GZIP', header=hdict, clobber=True)

    # return bitmask
    return None


bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/survey-bricks-dr9-{}.fits.gz'.format(field, field)))
# print(len(bricks))

# ########################## single brick ##########################
# mask = bricks['brickname']=='1092p320'
# brick_index = np.where(mask)[0][0]
# bitmask = get_pixel_bitmask(brick_index)
# ##################################################################

###################################################
if debug:
    np.random.seed(213)
    idx = np.random.choice(len(bricks), size=(10000), replace=False)
    bricks = bricks[idx]
###################################################

# random shuffle
np.random.seed(213)
bricks_list = np.random.choice(len(bricks), size=len(bricks), replace=False)
# split among the Cori nodes
bricks_list_split = np.array_split(bricks_list, n_task)
bricks_list = bricks_list_split[task_id]
print('Number of bricks in this node:', len(bricks_list))

time_start = time.time()
with Pool(processes=n_processes) as pool:
    res = pool.map(get_combined_bitmask, bricks_list, chunksize=1)

print('combine_lrg_masks {} {} {} Done!'.format(field, n_task, task_id), time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))

