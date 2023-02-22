from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from astropy.coordinates import SkyCoord
import healpy as hp


fn = '/global/cfs/cdirs/desi/users/rongpu/useful/lenz_hi/ebv_lhd.hpx.fits'
ebv = Table(fitsio.read(fn))['EBV']

ebv = hp.ud_grade(ebv, 2048, order_in='RING')

nside = 1024
pix_list = np.arange(hp.nside2npix(nside))
ra, dec = hp.pix2ang(nside, pix_list, nest=False, lonlat=True)
c = SkyCoord(ra, dec, unit='deg').galactic
l, b = c.l.to_value('deg'), c.b.to_value('deg')

ebv1 = hp.get_interp_val(ebv, l, b, nest=False, lonlat=True)
cat = Table()
cat['EBV'] = ebv1
cat.write('/global/cfs/cdirs/desi/users/rongpu/useful/lenz_hi/ebv_lhd_equatorial.hpx.fits')
