from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

custom_mask_fn = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/desi_custom_mask_v1.txt'

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

print('Custom circular mask:', len(circ_mask_data))
print('Custom rectangular mask:', len(rect_mask_data))

##############################################################

for field in ['north', 'south']:
    bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/survey-bricks-dr9-{}.fits.gz'.format(field, field)))

    mask_bricks = np.full(len(bricks), False)
    for radec in rect_mask_data:
        ramin, ramax, decmin, decmax = radec
        ramin -= 1.0
        decmin -= 0.5
        ramax += 1.0
        decmax += 0.5

        mask_bricks |= (bricks['ra']>ramin) & (bricks['ra']<ramax) & (bricks['dec']>decmin) & (bricks['dec']<decmax)

    print(field, np.sum(mask_bricks))
