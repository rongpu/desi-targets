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
custom_mask_fn = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/desi_custom_mask_v1.txt'

gaia_input_dir = '/global/cscratch1/sd/rongpu/desi/lrg_pixel_bitmask/dev1'
wise_input_dir = '/global/cscratch1/sd/rongpu/desi/lrg_pixel_bitmask/dev'
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

n_processes = 8


################################
debug = False
################################


def get_combined_bitmask(brick_index):

    brickname = str(bricks['brickname'][brick_index])

    output_path = os.path.join(output_dir, '{}/coadd/{}/{}/{}-lrgmask.fits.gz'.format(field, brickname[:3], brickname, brickname))
    # if os.path.isfile(output_path):
    #     return None

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

    gaiam_path = os.path.join(gaia_input_dir, '{}/coadd/{}/{}/{}-gaiamask.fits.gz'.format(field, brickname[:3], brickname, brickname))
    wisem_path = os.path.join(wise_input_dir, '{}/coadd/{}/{}/{}-wisemask.fits.gz'.format(field, brickname[:3], brickname, brickname))

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



####################### Add custom masks ########################

with open(custom_mask_fn, 'r') as f:
    lines = list(map(str.strip, f.readlines()))

# circular mask
ra, dec, radius = [], [], []
# rectangular mask
ramin, ramax, decmin, decmax = [], [], [], []

circ_mask_data = []
rect_mask_data = []

for line in lines:
    if line!='' and line[0]!='#':
        line = line[:line.find('#')]
        line = list(map(float, line.split(',')))
        if len(line)==3:
            circ_mask_data.append(line)
        elif len(line)==4:
            rect_mask_data.append(line)
        else:
            raise ValueError

circ_mask_data = np.array(circ_mask_data)
rect_mask_data = np.array(rect_mask_data)
# print('Custom circular mask:', len(circ_mask_data))
# print('Custom rectangular mask:', len(rect_mask_data))

##############################################################

bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/survey-bricks-dr9-{}.fits.gz'.format(field, field)))
# print(len(bricks))


mask_bricks = np.full(len(bricks), False)
for radec in rect_mask_data:
    ramin, ramax, decmin, decmax = radec
    ramin -= 1.0
    decmin -= 0.5
    ramax += 1.0
    decmax += 0.5

    mask_bricks |= (bricks['ra']>ramin) & (bricks['ra']<ramax) & (bricks['dec']>decmin) & (bricks['dec']<decmax)

print(field, np.sum(mask_bricks))
bricks = bricks[mask_bricks]

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

