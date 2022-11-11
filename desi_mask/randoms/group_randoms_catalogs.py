# Group randoms into fewer catalogs for more efficient queries

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

fns = glob.glob('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-[0-9]*')
fns.sort()
fns = [os.path.basename(tmp) for tmp in fns]
print(len(fns))

fns_split, file_indices_split = np.array_split(fns, 4), np.array_split(np.arange(len(fns)), 4)

for ii, fns in enumerate(fns_split):
    file_indices = file_indices_split[ii]
    cat_stack = []
    for jj, fn in enumerate(fns):
        print(ii, jj)
        cat = Table(fitsio.read('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'+fn, columns=['RA', 'DEC']))
        cat['file_index'] = np.array(file_indices[jj], dtype='int16')
        cat_stack.append(cat)
    cat_stack = vstack(cat_stack)
    cat_stack.write('/pscratch/sd/r/rongpu/tmp/randoms_group_{}.fits'.format(ii))
    del cat, cat_stack
