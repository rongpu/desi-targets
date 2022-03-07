# Create full sky EBV healpix map
# Run on NERSC

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
import healpy as hp

sys.path.append(os.path.expanduser('~/git/Python/useful'))
from get_ebv_from_map import get_ebv_from_map


nside_highres = 2048
npix = hp.nside2npix(nside_highres)

ra, dec = hp.pix2ang(nside_highres, np.arange(npix), lonlat=True, nest=False)
ebv_highres = get_ebv_from_map([ra, dec])

for nside in [32, 64, 128, 256, 512, 1024]:
    print(nside)
    ebv = hp.pixelfunc.ud_grade(ebv_highres, nside, order_in='RING', order_out='RING')
    np.save('/global/cfs/cdirs/desi/users/rongpu/useful/sfd_healpix/sfd_ebv_{}_ring.npy'.format(nside), ebv)


