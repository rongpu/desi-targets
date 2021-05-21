# Minimalist GAIA catalog for testing DESI selection around bright stars

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

import healpy as hp
from astropy import units as u
from astropy.coordinates import SkyCoord

galactic_b_limit_pos, galactic_b_limit_neg = 14., -20.

gaia_dir = '/project/projectdirs/cosmo/data/gaia/dr2/healpix'
output_path = '/global/cfs/cdirs/desi/users/rongpu/useful/gaia_dr2_g_18.fits'

# gaia_dir = '/project/projectdirs/cosmo/data/gaia/edr3/healpix'
# output_path = '/global/cfs/cdirs/desi/users/rongpu/useful/gaia_edr3_g_18.fits'

gaia_nside = 32
print('Healpix pixel size (square deg): {:.5f}'.format(hp.nside2pixarea(gaia_nside, degrees=True)))

gaia_npix = hp.nside2npix(gaia_nside)
gaia_hp_ra, gaia_hp_dec = hp.pix2ang(gaia_nside, np.arange(gaia_npix), nest=True, lonlat=True)
c = SkyCoord(gaia_hp_ra, gaia_hp_dec, unit='deg').galactic
gaia_hp_l, gaia_hp_b = c.l.to_value('deg'), c.b.to_value('deg')

mask = (gaia_hp_b>galactic_b_limit_pos-2) | (gaia_hp_b<galactic_b_limit_neg+2)
mask &= (gaia_hp_dec>-32)
gaia_list = np.where(mask)[0]
print(len(gaia_list))

gaia = []
for index, hp_index in enumerate(gaia_list):
    gaia_fn = str(hp_index).zill(5)
    if index%100==0:
        print('{}/{}, {}'.format(index, len(gaia_list), gaia_fn))
    tmp = Table(fitsio.read(os.path.join(gaia_dir, 'healpix-{}.fits'.format(gaia_fn)), columns=['SOURCE_ID', 'RA', 'DEC', 'PHOT_G_MEAN_MAG', 'PHOT_BP_MEAN_MAG', 'PHOT_RP_MEAN_MAG', 'PHOT_G_MEAN_FLUX_OVER_ERROR']))
    mask = tmp['PHOT_G_MEAN_MAG']<18.0
    mask &= tmp['PHOT_G_MEAN_FLUX_OVER_ERROR']>0
    mask &= (tmp['DEC']>-32)
    if np.sum(mask)>0:
        c = SkyCoord(tmp['RA'], tmp['DEC'], unit='deg').galactic
        gaia_hp_l, gaia_hp_b = c.l.to_value('deg'), c.b.to_value('deg')
        mask &= (gaia_hp_b>galactic_b_limit_pos) | (gaia_hp_b<galactic_b_limit_neg)
        if np.sum(mask)>0:
            tmp = tmp[mask]
            gaia.append(tmp)
gaia = vstack(gaia)
print(len(gaia))

gaia.write(output_path)
