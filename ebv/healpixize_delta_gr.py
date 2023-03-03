from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

from multiprocessing import Pool
import healpy as hp

nside = 64
n_processes = 128


def do_stuff(pix_idx):

    pix_list = pix_unique[pix_idx]

    hp_table = Table()
    hp_table['HPXPIXEL'] = pix_list
    # hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_list, nest=False, lonlat=True)

    hp_table['delta_gr_mean'] = 0.
    hp_table['delta_gr_median'] = 0.
    hp_table['n_star'] = 0
    hp_table['EBV'] = 0.

    for index in np.arange(len(pix_idx)):
        idx = pixorder[pixcnts[pix_idx[index]]:pixcnts[pix_idx[index]+1]]
        v = cat['DATA_G-R'][idx]-cat['MODEL_G-R'][idx]
        hp_table['delta_gr_mean'][index] = np.mean(v)
        hp_table['delta_gr_median'][index] = np.median(v)
        hp_table['n_star'][index] = len(v)
        hp_table['EBV'][index] = np.mean(cat['EBV'][idx])

    return hp_table


cat = Table(fitsio.read('/pscratch/sd/r/rongpu/ebv/desi_std/desi_standard_stars_iron.fits'))

mask = cat['FIBERSTATUS']==0
cat = cat[mask]
print(len(cat))

# mask = cat['BLUE_SNR']>10
# cat = cat[mask]
# print(len(cat))

exp = Table(fitsio.read('/global/homes/j/jguy/redux/iron/exposures-iron.fits', ext='EXPOSURES'))
mask = (exp['SURVEY']=='sv1') | (exp['SURVEY']=='sv3') | (exp['SURVEY']=='main')
mask = np.in1d(cat['EXPID'], exp['EXPID'][mask])
cat = cat[mask]
print(len(cat))

mask = np.isfinite(cat['DATA_G-R']-cat['MODEL_G-R'])
cat = cat[mask]
print(len(cat))

npix = hp.nside2npix(nside)

pix_allobj = hp.pixelfunc.ang2pix(nside, cat['TARGET_RA'], cat['TARGET_DEC'], lonlat=True)
pix_unique, pixcnts = np.unique(pix_allobj, return_counts=True)

pixcnts = np.insert(pixcnts, 0, 0)
pixcnts = np.cumsum(pixcnts)

pixorder = np.argsort(pix_allobj)

# split among the Cori processors
pix_idx_split = np.array_split(np.arange(len(pix_unique)), n_processes)

# start multiple worker processes
with Pool(processes=n_processes) as pool:
    res = pool.map(do_stuff, pix_idx_split)

hp_table = vstack(res)
hp_table.sort('HPXPIXEL')

hp_table.write('/pscratch/sd/r/rongpu/ebv/desi_std/delta_gr_sv1sv3main_nside_{}.fits'.format(nside), overwrite=True)
