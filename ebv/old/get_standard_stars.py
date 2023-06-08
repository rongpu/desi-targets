from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from multiprocessing import Pool

exp = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/iron/exposures-iron.fits', ext='EXPOSURES'))
print(len(exp))


def get_std(exp_index):

    cat_stack = []

    for petal_loc in range(10):

        night = exp['NIGHT'][exp_index]
        expid = exp['EXPID'][exp_index]
        expid_str = str(expid).zfill(8)
        fn = f'/global/cfs/cdirs/desi/spectro/redux/iron/exposures/{night}/{expid_str}/stdstars-{petal_loc}-{expid_str}.fits.gz'

        if not os.path.isfile(fn):
            continue

        cat1 = Table(fitsio.read(fn, ext='METADATA'))
        cat2 = Table(fitsio.read(fn, ext='FIBERMAP'))
        # cat1.sort('TARGETID')
        # cat2.sort('TARGETID')
        assert np.all(cat1['TARGETID']==cat2['TARGETID']) and np.all(cat1['FIBER']==cat2['FIBER'])

        cat2.remove_columns(['TARGETID', 'FIBER'])
        cat = hstack([cat1, cat2])
        # print(expid, len(cat))

        cat_stack.append(cat)

    if len(cat_stack)==0:
        return None
    else:
        cat_stack = vstack(cat_stack)
        cat_stack['EXPID'] = expid
        return cat_stack


n_processess = 128
with Pool(processes=n_processess) as pool:
    res = pool.map(get_std, np.arange(len(exp)))

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

cat = vstack(res)
cat.write('/pscratch/sd/r/rongpu/ebv/desi_std/desi_standard_stars_iron.fits')
