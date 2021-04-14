from __future__ import division, print_function
import numpy as np

# Anand Raichoor's code for checking if ra, dec are in the DES footprint
def get_isdes(ra, dec, nside):
    import pymangle
    import healpy as hp
    
    npix = hp.nside2npix(nside)
    # checking hp pixels
    mng = pymangle.Mangle('des.ply')
    theta, phi = hp.pix2ang(nside, np.arange(npix), nest=False)
    hpra, hpdec= 180./np.pi*phi, 90.-180./np.pi*theta
    hpindes = (mng.polyid(hpra, hpdec)!=-1).astype(int)
    # pixels with all neighbours in des
    hpindes_secure = np.array([i for i in range(npix) if hpindes[i]+hpindes[hp.get_all_neighbours(nside, i)].sum()==9])
    # pixels with all neighbours outside des
    hpoutdes_secure = np.array([i for i in range(npix) if hpindes[i]+hpindes[hp.get_all_neighbours(nside, i)].sum()==0])
    # hpind to be checked
    tmp = np.ones(npix, dtype=bool)
    tmp[hpindes_secure] = False
    tmp[hpoutdes_secure]= False
    hp_tbc = np.arange(npix)[tmp]

    # now checking indiv. obj. in the tbc pixels
    hppix = hp.ang2pix(nside, (90.-dec)*np.pi/180., ra*np.pi/180., nest=False)
    hpind = np.unique(hppix)

    isdes = np.zeros(len(ra), dtype=bool)
    isdes[np.in1d(hppix, hpindes_secure)] = True
    tbc = np.where(np.in1d(hppix, hp_tbc))[0]
    tbcisdes = (mng.polyid(ra[tbc], dec[tbc])!=-1)
    isdes[tbc][tbcisdes] = True

    return isdes