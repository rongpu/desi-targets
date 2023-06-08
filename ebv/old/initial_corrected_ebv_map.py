from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

for nside in [32, 64, 128, 256]:
    maps = Table(fitsio.read('/pscratch/sd/r/rongpu/ebv/desi_std/delta_gr_sv1sv3main_nside_{}.fits'.format(nside)))
    maps.rename_column('EBV', 'EBV_SFD')
    maps['EBV_NEW'] = maps['EBV_SFD'] + 1/(3.214-2.165) * maps['delta_gr_mean'] - 0.023
    maps.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/test/initial_corrected_ebv_map_nside_{}.fits'.format(nside), overwrite=True)
