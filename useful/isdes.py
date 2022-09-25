from __future__ import division, print_function
import os
import numpy as np
from astropy.table import Table, vstack, hstack
import healpy as hp


# # Anand Raichoor's code for checking if ra, dec are in the DES footprint
# import pymangle
# for nside in [32, 64, 128, 256, 512, 1024]:
#     npix = hp.nside2npix(nside)
#     hp_ra, hp_dec = hp.pix2ang(nside, np.arange(npix), nest=False, lonlat=True)

#     # checking hp pixels
#     mng = pymangle.Mangle(os.path.expanduser('~/git/desi-targets/useful/in_des/des.ply'))
#     theta, phi = hp.pix2ang(nside, np.arange(npix), nest=False)
#     hp_ra, hpd_ec= 180./np.pi*phi, 90.-180./np.pi*theta
#     indes = mng.polyid(hp_ra, hpd_ec)!=-1

#     cat = Table()
#     cat['in_des'] = indes
#     cat.write(os.path.expanduser('~/git/desi-targets/useful/in_des/hp_in_des_{}_ring.fits.gz'.format(nside)))


def get_isdes(ra, dec, nside=256):

    hpix = hp.ang2pix(nside, ra, dec, nest=False, lonlat=True)
    cat = Table.read(os.path.expanduser('~/git/desi-targets/useful/in_des/hp_in_des_{}_ring.fits.gz'.format(nside)))
    isdes = np.array(cat['in_des'][hpix])

    return isdes

