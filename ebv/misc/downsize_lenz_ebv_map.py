from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from astropy.coordinates import SkyCoord
import healpy as hp

sys.path.append(os.path.expanduser('~/git/desi-examples/misc/misc'))
from healpix_util import downsize_hp_map


nside_in = 1024
fn = '/global/cfs/cdirs/desi/users/rongpu/useful/lenz_hi/ebv_lhd.hpx.fits'
ebv = Table(fitsio.read(fn))
ebv['HPXPIXEL'] = np.arange(hp.nside2npix(nside_in))

stats_dict = {'EBV': np.nanmean}
columns = ['HPXPIXEL', 'EBV']
for nside_out in [64, 128, 256, 512]:
    ebv_new = downsize_hp_map(1024, nside_out, ebv[columns], stats_dict=stats_dict, weights=None, n_processes=128)
    ebv_new.write('/global/cfs/cdirs/desi/users/rongpu/useful/lenz_hi/ebv_lhd_{}.hpx.fits'.format(nside_out), overwrite=True)
