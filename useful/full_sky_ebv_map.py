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

from astropy.coordinates import SkyCoord

sys.path.append(os.path.expanduser('~/git/Python/useful'))
from get_ebv_from_map import get_ebv_from_map


nside_highres = 2048
npix = hp.nside2npix(nside_highres)

# Equatorial coordinates
ra, dec = hp.pix2ang(nside_highres, np.arange(npix), lonlat=True, nest=False)
ebv_highres = get_ebv_from_map([ra, dec])
np.save('/global/cfs/cdirs/desi/users/rongpu/useful/sfd_healpix/sfd_ebv_{}_ring.npy'.format(nside_highres), ebv_highres)
for nside in [32, 64, 128, 256, 512, 1024]:
    print(nside)
    ebv = hp.pixelfunc.ud_grade(ebv_highres, nside, order_in='RING', order_out='RING')
    output_path = '/global/cfs/cdirs/desi/users/rongpu/useful/sfd_healpix/sfd_ebv_{}_ring.npy'.format(nside)
    if not os.path.isfile(output_path):
        np.save(output_path, ebv)

# Galactic coordinates
l, b = hp.pix2ang(nside_highres, np.arange(npix), lonlat=True, nest=False)
c = SkyCoord(l, b, unit='deg', frame='galactic').icrs
ra, dec = c.ra.to_value('deg'), c.dec.to_value('deg')
ebv_highres = get_ebv_from_map([ra, dec])
np.save('/global/cfs/cdirs/desi/users/rongpu/useful/sfd_healpix/sfd_ebv_{}_galactic_ring.npy'.format(nside_highres), ebv_highres)
for nside in [32, 64, 128, 256, 512, 1024]:
    print(nside)
    ebv = hp.pixelfunc.ud_grade(ebv_highres, nside, order_in='RING', order_out='RING')
    output_path = '/global/cfs/cdirs/desi/users/rongpu/useful/sfd_healpix/sfd_ebv_{}_galactic_ring.npy'.format(nside)
    if not os.path.isfile(output_path):
        np.save(output_path, ebv)



