# Create symlinks for sweep catalogs -- use ZP-corrected catalogs when available

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

sweep_dir_orig = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0'
sweep_dir_corr = '/pscratch/sd/r/rongpu/dr9_desi_photoz/sweep_zp_corrected'
symlink_dir = '/global/cfs/cdirs/desicollab/users/rongpu/data/dr9/sweep_symlinks'

sweep_all_path = sorted(glob.glob(os.path.join(sweep_dir_orig, '*.fits')))
sweep_fn_list = [os.path.basename(sweep_all_path[ii]) for ii in range(len(sweep_all_path))]

for sweep_fn in sweep_fn_list:

    symlink_path = os.path.join(symlink_dir, sweep_fn)

    sweep_path = os.path.join(sweep_dir_corr, sweep_fn)
    if not os.path.isfile(sweep_path):
        sweep_path = os.path.join(sweep_dir_orig, sweep_fn)

    os.symlink(sweep_path, symlink_path)
    print(symlink_path)
