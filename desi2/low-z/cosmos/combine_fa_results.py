from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

fns = sorted(glob.glob('/global/cfs/cdirs/desi/users/raichoor/fiberassign-cosmoslowz/20220404/999/*gz'))

cat = []
for fn in fns:
    tmp = Table(fitsio.read(fn, ext="FIBERASSIGN"))
    tmp['tileid'] = int(fn[-14:-8])
    cat.append(tmp)
cat = vstack(cat)
cat.write("/global/cfs/cdirs/desi/users/rongpu/tmp/tmp-cosmoslowz-20220404-fiberassign.fits", overwrite=False)
