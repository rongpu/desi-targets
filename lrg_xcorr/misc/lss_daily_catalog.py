# Trim and clean up the LSS daily catalog

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack
import fitsio

# cat = Table(fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/daily/LSScats/test/LRG_full.dat.fits'))
# cat.rename_column('Z_not4clus', 'Z')

cat = Table(fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/daily/LSScats/test/LRG_full_noveto.dat.fits'))

print(len(cat), len(np.unique(cat['TARGETID'])))

# cat['EFFTIME_LRG'] = 12.15 * cat['TSNR2_LRG']

# mask = cat['COADD_FIBERSTATUS']==0
# cat = cat[mask]
# print(len(cat), len(np.unique(cat['TARGETID'])))

# mask = cat['EFFTIME_LRG']>800
# print(np.sum(mask)/len(mask))
# cat = cat[mask]
# print(len(cat), len(np.unique(cat['TARGETID'])))

# # Redshift quality cut
# cat['q'] = cat['ZWARN']==0
# cat['q'] &= cat['Z']<1.5
# cat['q'] &= cat['DELTACHI2']>15
# print(np.sum(cat['q']), np.sum(cat['q'])/len(cat))

# cat['isstar'] = (cat['SPECTYPE']=='STAR') | (cat['Z']<0.0003)
# print(np.sum(cat['isstar'])/len(cat))

# columns_to_keep = ['TARGETID', 'SUBPRIORITY', 'PRIORITY_INIT', 'TARGET_STATE', 'TILEID', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'DELTACHI2', 'FIBER', 'COADD_FIBERSTATUS', 'FIBERASSIGN_X', 'FIBERASSIGN_Y', 'PRIORITY', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'MEAN_DELTA_X', 'RMS_DELTA_X', 'MEAN_DELTA_Y', 'RMS_DELTA_Y', 'MEAN_PSF_TO_FIBER_SPECFLUX', 'TSNR2_ELG', 'TSNR2_LYA', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG', 'ZWARN_MTL', 'GOODHARDLOC', 'NTILE', 'GOODTSNR', 'RA', 'DEC', 'EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'MASKBITS', 'PHOTSYS', 'DESI_TARGET', 'BGS_TARGET', 'COMP_TILE', 'FRACZ_TILELOCID', 'lrg_mask', 'q', 'isstar']
columns_to_keep = ['TARGETID', 'SUBPRIORITY', 'PRIORITY_INIT', 'TARGET_STATE', 'TILEID', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'DELTACHI2', 'FIBER', 'COADD_FIBERSTATUS', 'FIBERASSIGN_X', 'FIBERASSIGN_Y', 'PRIORITY', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'MEAN_DELTA_X', 'RMS_DELTA_X', 'MEAN_DELTA_Y', 'RMS_DELTA_Y', 'MEAN_PSF_TO_FIBER_SPECFLUX', 'TSNR2_ELG', 'TSNR2_LYA', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG', 'ZWARN_MTL', 'GOODHARDLOC', 'NTILE', 'GOODTSNR', 'RA', 'DEC', 'EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'MASKBITS', 'PHOTSYS', 'DESI_TARGET', 'BGS_TARGET', 'COMP_TILE', 'FRACZ_TILELOCID', 'lrg_mask']
cat = cat[columns_to_keep]

# cat.write('/global/cfs/cdirs/desi/users/rongpu/data/ls_dr9.0_photoz/misc/lss_daily_y1_lrg_full_20221031.fits', overwrite=True)
cat.write('/global/cfs/cdirs/desi/users/rongpu/data/ls_dr9.0_photoz/misc/lss_daily_y1_lrg_full_noveto_20221031.fits', overwrite=True)

