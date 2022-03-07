from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

from desimodel.footprint import is_point_in_desi
from astropy.table import Table


def get_isdesi(tilesfn, ra, dec):
    tiles = Table.read(tilesfn)
    tiles = tiles[(tiles["PROGRAM"] == "DARK") & (tiles["IN_DESI"])]
    ii = is_point_in_desi(tiles, ra, dec)
    isdesi = np.zeros(len(ra), dtype=bool)
    isdesi[ii] = True
    return isdesi


tilesfn = "/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-main.ecsv"

data = {}

for nside in [128, 256, 512, 1024]:

    npix = hp.nside2npix(nside)
    hp_idx = np.arange(npix)
    ra, dec = hp.pix2ang(nside, hp_idx, lonlat=True)

    indesi = get_isdesi(tilesfn, ra, dec)

    data[str(nside)] = hp_idx[indesi].copy()

np.save('/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/misc/in_desi_healpix_ring.npy', data)

# hp_in_desi = np.load('/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/misc/in_desi_healpix_ring.npy', allow_pickle=True).item()[str(nside)]
