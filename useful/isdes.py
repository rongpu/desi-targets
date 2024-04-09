from __future__ import division, print_function
import os
import numpy as np
from astropy.table import Table, vstack, hstack
import healpy as hp
import fitsio


# # Anand Raichoor's code for checking if ra, dec are in the DES footprint
# import pymangle
# for nside in [32, 64, 128, 256, 512, 1024]:
#     npix = hp.nside2npix(nside)
#     hp_ra, hp_dec = hp.pix2ang(nside, np.arange(npix), nest=False, lonlat=True)

#     # checking hp pixels
#     mng = pymangle.Mangle('/global/cfs/cdirs/desi/users/rongpu/useful/in_des/des.ply')
#     theta, phi = hp.pix2ang(nside, np.arange(npix), nest=False)
#     hp_ra, hpd_ec= 180./np.pi*phi, 90.-180./np.pi*theta
#     indes = mng.polyid(hp_ra, hpd_ec)!=-1

#     cat = Table()
#     cat['in_des'] = indes
#     cat.write('/global/cfs/cdirs/desi/users/rongpu/useful/in_des/hp_in_des_{}_ring.fits.gz'.format(nside))

# for nside in [32, 64, 128, 256, 512, 1024]:
#     cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/useful/in_des/hp_in_des_{}_ring.fits.gz'.format(nside)))
#     cat['HPXPIXEL'] = np.arange(len(cat))
#     cat['HPXPIXEL'] = hp.ring2nest(nside, cat['HPXPIXEL'])
#     cat.sort('HPXPIXEL')
#     cat.remove_column('HPIXPIXEL')
#     cat.write('/global/cfs/cdirs/desi/users/rongpu/useful/in_des/hp_in_des_{}_nest.fits.gz'.format(nside))


def get_isdes(ra, dec, nside=256):

    hpix = hp.ang2pix(nside, ra, dec, nest=False, lonlat=True)
    cat = Table.read('/global/cfs/cdirs/desi/users/rongpu/useful/in_des/hp_in_des_{}_ring.fits.gz'.format(nside))
    isdes = np.array(cat['in_des'][hpix])

    return isdes

