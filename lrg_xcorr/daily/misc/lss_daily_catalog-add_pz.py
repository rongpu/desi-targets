from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack
import fitsio

# cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/ls_dr9.0_photoz/misc/lss_daily_y1_lrg_full_20221031.fits'))
cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/ls_dr9.0_photoz/misc/lss_daily_y1_lrg_full_noveto_20221031.fits'))

sys.path.append(os.path.expanduser('~/git/desi-examples/misc'))
from misc.get_sweep_columns import get_sweep_columns

columns = ['Z_PHOT_MEAN', 'Z_PHOT_MEDIAN', 'Z_PHOT_STD', 'Z_PHOT_L68', 'Z_PHOT_U68', 'Z_PHOT_L95', 'Z_PHOT_U95', 'Z_SPEC', 'SURVEY', 'TRAINING']

pz_dir = '/global/cfs/cdirs/desi/users/rongpu/data/ls_dr9.0_desi_photoz/pz/'
pz = get_sweep_columns(cat, columns, n_processes=128, pz_dir=pz_dir)
pz.write('/global/cfs/cdirs/desi/users/rongpu/data/ls_dr9.0_photoz/misc/lss_daily_y1_lrg_full_20221031-pz_new.fits')
# pz.write('/global/cfs/cdirs/desi/users/rongpu/data/ls_dr9.0_photoz/misc/lss_daily_y1_lrg_full_noveto_20221031-pz_new.fits')

# pz = get_sweep_columns(cat, columns, n_processes=128, pz_dir=None)
# pz.write('/global/cfs/cdirs/desi/users/rongpu/data/ls_dr9.0_photoz/misc/lss_daily_y1_lrg_full_20221031-pz_dr9.fits')

