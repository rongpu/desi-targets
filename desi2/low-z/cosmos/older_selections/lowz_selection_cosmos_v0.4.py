from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

# COSMOS field center
ra, dec = 150.11917, 2.20583
ramin, ramax, decmin, decmax = ra - 2, ra + 2, dec - 2, dec + 2

# ############################ Start with sweep catalogs ############################

# sweep_fns = ['sweep-140p000-150p005.fits', 'sweep-150p000-160p005.fits']
# print(len(sweep_fns))

# columns = ['RELEASE', 'BRICKID', 'OBJID', 'RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'TYPE', 'EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG', 'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z']

# cat_stack = []
# for sweep_fn in sweep_fns:
#     cat = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/'+sweep_fn, columns=columns))
#     with warnings.catch_warnings():
#         warnings.simplefilter("ignore")
#         zmag = 22.5 - 2.5*np.log10(cat['FLUX_Z']) - 1.211 * cat['EBV']
#         zfibermag = 22.5 - 2.5*np.log10(cat['FIBERFLUX_Z']) - 1.211 * cat['EBV']
#     mask = (zmag<21.5) | (zfibermag<22.5)
#     mask &= np.isfinite(zfibermag)
#     mask &= np.isfinite(zmag)
#     mask &= (cat['NOBS_G']>=1) & (cat['NOBS_R']>=1) & (cat['NOBS_Z']>=1)
#     print(sweep_fn, np.sum(mask)/len(mask), np.sum(mask))
#     mask &= (cat['RA']>ramin) & (cat['RA']<ramax) & (cat['DEC']>decmin) & (cat['DEC']<decmax)
#     print(sweep_fn, np.sum(mask)/len(mask), np.sum(mask))
#     cat = cat[mask]
#     cat_stack.append(cat)
# cat = vstack(cat_stack)
# print(len(cat))

# cat.write('/global/cscratch1/sd/rongpu/temp/sweep_cosmos_secondary.fits', overwrite=True)

############################ Basic quality cuts ############################

cat = Table.read('/global/cscratch1/sd/rongpu/temp/sweep_cosmos_secondary.fits')
print(len(cat))

from desitarget.targets import decode_targetid, encode_targetid
cat['TARGETID'] = encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE'])

min_nobs = 1
mask = (cat['NOBS_G']>=min_nobs) & (cat['NOBS_R']>=min_nobs) & (cat['NOBS_Z']>=min_nobs)
cat = cat[mask]
print(len(cat))

maskbits = [1, 12, 13]
mask = np.full(len(cat), True)
for bit in maskbits:
    mask &= (cat['MASKBITS'] & 2**bit)==0
cat = cat[mask]
print(len(cat))

mask_quality = np.full(len(cat), True)
mask_quality &= (cat['FLUX_R'] > 0)   # quality in r.
mask_quality &= (cat['FLUX_Z'] > 0) & (cat['FIBERFLUX_Z'] > 0)   # quality in z.
mask_quality &= (cat['FLUX_W1'] > 0)  # quality in W1.
mask_quality &= (cat['GAIA_PHOT_G_MEAN_MAG'] == 0) | (cat['GAIA_PHOT_G_MEAN_MAG'] > 18)  # remove bright GAIA sources
mask_quality &= cat['FIBERTOTFLUX_R'] < 10**(-0.4*(18.0-22.5))  # remove objects with rfibertot < 18.0
mask_quality &= cat['FIBERTOTFLUX_Z'] < 10**(-0.4*(17.5-22.5))  # remove objects with zfibertot < 17.5
cat = cat[mask_quality]
print(len(cat))

from astropy.coordinates import SkyCoord
from astropy import units as u

search_radius = 2.0
sky1 = SkyCoord(ra*u.degree, dec*u.degree, frame='icrs')
sky2 = SkyCoord(cat['RA']*u.degree, cat['DEC']*u.degree, frame='icrs')
mask = np.array(sky2.separation(sky1).to(u.degree))<search_radius
cat = cat[mask]
print(len(cat))

stats_search_radius = 1.62
sky1 = SkyCoord(ra*u.degree, dec*u.degree, frame='icrs')
sky2 = SkyCoord(cat['RA']*u.degree, cat['DEC']*u.degree, frame='icrs')
mask = np.array(sky2.separation(sky1).to(u.degree))<stats_search_radius
cat['in_cosmos'] = mask.copy()
print(np.sum(cat['in_cosmos']))

############################ Selection cuts ############################

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    cat['gmag'] = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_G']*10**(0.4*3.214*cat['EBV']), 1e-7, None))
    cat['rmag'] = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_R']*10**(0.4*2.165*cat['EBV']), 1e-7, None))
    cat['zmag'] = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_Z']*10**(0.4*1.211*cat['EBV']), 1e-7, None))
    cat['w1mag'] = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_W1']*10**(0.4*0.184*cat['EBV']), 1e-7, None))
    cat['w2mag'] = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_W2']*10**(0.4*0.113*cat['EBV']), 1e-7, None))
    cat['gfibermag'] = 22.5 - 2.5*np.log10(np.clip(cat['FIBERFLUX_G']*10**(0.4*3.214*cat['EBV']), 1e-7, None))
    cat['rfibermag'] = 22.5 - 2.5*np.log10(np.clip(cat['FIBERFLUX_R']*10**(0.4*2.165*cat['EBV']), 1e-7, None))
    cat['zfibermag'] = 22.5 - 2.5*np.log10(np.clip(cat['FIBERFLUX_Z']*10**(0.4*1.211*cat['EBV']), 1e-7, None))

primary = cat['zfibermag']<21.6
primary &= cat['zmag'] - cat['w1mag'] > 0.8 * (cat['rmag'] - cat['zmag']) - 1.1  # LRG non-stellar cut
print(np.sum(primary), np.sum(primary)/len(cat))

filler = (cat['zfibermag']>=21.6) & (cat['zfibermag']<22.4)
print(np.sum(filler), np.sum(filler)/len(cat))

cat['primary'] = primary.copy()
cat['filler'] = filler.copy()

mask = cat['primary'] | cat['filler']
cat = cat[mask]
print(len(cat))

#################### Flag existing DESI redshifts ######################

# t1 = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv1-dark.fits'))
# mask = (t1['TARGET_RA']>ramin) & (t1['TARGET_RA']<ramax) & (t1['TARGET_DEC']>decmin) & (t1['TARGET_DEC']<decmax)
# t1 = t1[mask]
# t1['bright'] = False
# t2 = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv3-dark.fits'))
# mask = (t2['TARGET_RA']>ramin) & (t2['TARGET_RA']<ramax) & (t2['TARGET_DEC']>decmin) & (t2['TARGET_DEC']<decmax)
# t2 = t2[mask]
# t2['bright'] = False

# t3 = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv1-bright.fits'))
# mask = (t3['TARGET_RA']>ramin) & (t3['TARGET_RA']<ramax) & (t3['TARGET_DEC']>decmin) & (t3['TARGET_DEC']<decmax)
# t3 = t3[mask]
# t3['bright'] = True
# t4 = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv3-bright.fits'))
# mask = (t4['TARGET_RA']>ramin) & (t4['TARGET_RA']<ramax) & (t4['TARGET_DEC']>decmin) & (t4['TARGET_DEC']<decmax)
# t4 = t4[mask]
# t4['bright'] = True

# # t5 = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/guadalupe/zcatalog/zpix-main-bright.fits'))
# # mask = (t5['TARGET_RA']>ramin) & (t5['TARGET_RA']<ramax) & (t5['TARGET_DEC']>decmin) & (t5['TARGET_DEC']<decmax)
# # t5 = t5[mask]
# # t5['bright'] = False

# # obs = vstack([t1, t2, t3, t3, t4, t5], join_type='inner')
# obs = vstack([t1, t2, t3, t4], join_type='inner')
# print(len(obs), len(np.unique(obs['TARGETID'])))

# obs['EFFTIME_LRG'] = 12.15 * obs['TSNR2_LRG']

# # Remove duplicates keeping the highest EFFTIME_LRG
# obs.sort('EFFTIME_LRG', reverse=True)
# _, idx_keep = np.unique(obs['TARGETID'], return_index=True)
# obs = obs[idx_keep]
# print(len(obs), len(np.unique(obs['TARGETID'])))

# obs.rename_columns(['TARGET_RA', 'TARGET_DEC'], ['RA', 'DEC'])
# obs.write('/global/cscratch1/sd/rongpu/temp/fuji_in_cosmos.fits')

obs = Table(fitsio.read('/global/cscratch1/sd/rongpu/temp/fuji_in_cosmos.fits'))

mask = obs['ZWARN']==0
mask &= obs['COADD_FIBERSTATUS']==0
# mask &= ((obs['bright']==False) & obs['DELTACHI2']>15) & ((obs['bright']==True) & (obs['DELTACHI2']>40))
mask &= (obs['DELTACHI2']>40)
obs = obs[mask]
print(len(obs))

mask_obs = np.in1d(cat['TARGETID'], obs['TARGETID'])
cat['observed'] = mask_obs.copy()

#################### Exclude SV3 targets ######################

# sv3_dir = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/0.57.0'

# sv3lrg = []
# for field in ['north', 'south']:
#     tmp = Table(fitsio.read(os.path.join(sv3_dir, 'dr9_sv3_lrg_{}_0.57.0_basic.fits'.format(field)), columns=['TARGETID']))
#     sv3lrg.append(hstack([tmp]))
# sv3lrg = vstack(sv3lrg)
# print(len(sv3lrg))

# sv3bgs = []
# for field in ['north', 'south']:
#     tmp = Table(fitsio.read(os.path.join(sv3_dir, 'dr9_sv3_bgs_any_{}_0.57.0_basic.fits'.format(field)), columns=['TARGETID', 'SV3_BGS_TARGET']))
#     sv3bgs.append(hstack([tmp]))
# sv3bgs = vstack(sv3bgs)
# print(len(sv3bgs))
# # Select BGS Bright and BGS Faint targets
# mask = (sv3bgs['SV3_BGS_TARGET'] & 2**1 > 0) | (sv3bgs['SV3_BGS_TARGET'] & 2**0 > 0)
# sv3bgs = sv3bgs[mask]
# print(len(sv3bgs))

# sv3qso = []
# for field in ['north', 'south']:
#     tmp = Table(fitsio.read(os.path.join(sv3_dir, 'dr9_sv3_qso_{}_0.57.0_basic.fits'.format(field)), columns=['TARGETID']))
#     sv3qso.append(hstack([tmp]))
# sv3qso = vstack(sv3qso)
# print(len(sv3qso))

# cat['sv3lrg'] = np.in1d(cat['TARGETID'], sv3lrg['TARGETID'])
# cat['sv3bgs'] = np.in1d(cat['TARGETID'], sv3bgs['TARGETID'])
# cat['sv3qso'] = np.in1d(cat['TARGETID'], sv3qso['TARGETID'])

import healpy as hp

nside = 8
pix_list = np.unique(hp.ang2pix(nside, cat['RA'], cat['DEC'], nest=True, lonlat=True))

bright = []
dark = []
for pix in pix_list:
    bright_path = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.57.0/targets/sv3/resolve/bright/sv3targets-bright-hp-{}.fits'.format(pix)
    dark_path = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.57.0/targets/sv3/resolve/dark/sv3targets-dark-hp-{}.fits'.format(pix)
    bright.append(Table(fitsio.read(bright_path, columns=['TARGETID', 'SV3_BGS_TARGET'])))
    dark.append(Table(fitsio.read(dark_path, columns=['TARGETID', 'SV3_DESI_TARGET'])))
bright = vstack(bright)
dark = vstack(dark)

# Select BGS Bright and BGS Faint targets
mask = (bright['SV3_BGS_TARGET'] & 2**1 > 0) | (bright['SV3_BGS_TARGET'] & 2**0 > 0)
bright = bright[mask]

# Select LRG and QSO targets
mask = (dark['SV3_DESI_TARGET'] & 2**0 > 0) | (dark['SV3_DESI_TARGET'] & 2**2 > 0)
dark = dark[mask]
dark['lrg'] = (dark['SV3_DESI_TARGET'] & 2**0 > 0)
dark['qso'] = (dark['SV3_DESI_TARGET'] & 2**2 > 0)

cat['sv3lrg'] = np.in1d(cat['TARGETID'], dark['TARGETID'][dark['lrg']])
cat['sv3qso'] = np.in1d(cat['TARGETID'], dark['TARGETID'][dark['qso']])
cat['sv3bgs'] = np.in1d(cat['TARGETID'], bright['TARGETID'])

cat['is_target'] = (~cat['sv3lrg']) & (~cat['sv3qso']) & (~(cat['sv3bgs'] & (cat['rfibermag']<21.)))

print(np.sum(cat['is_target'] & cat['in_cosmos']))

cat.remove_columns(['in_cosmos', 'gmag', 'rmag', 'zmag', 'w1mag', 'w2mag', 'rfibermag', 'zfibermag'])

# Assign priorities
cat['NUMOBS_INIT'] = 0
cat['NUMOBS_INIT'][cat['primary']] = 1
cat['NUMOBS_INIT'][cat['filler']] = 4
cat['PRIORITY_INIT'] = 0
cat['PRIORITY_INIT'][cat['primary']] = 9500
cat['PRIORITY_INIT'][cat['filler']] = 9000
cat['PRIORITY_DONE'] = 0
cat['PRIORITY_DONE'][cat['primary']] = 8500
cat['PRIORITY_DONE'][cat['filler']] = 8000
cat.write('/global/cfs/cdirs/desi/users/rongpu/misc/cosmos_lowz_sample_v0.4_scheme_1.fits')

# Assign priorities
cat['NUMOBS_INIT'] = 0
cat['NUMOBS_INIT'][cat['primary']] = 1
cat['NUMOBS_INIT'][cat['filler']] = 1
cat['PRIORITY_INIT'] = 0
cat['PRIORITY_INIT'][cat['primary']] = 9500
cat['PRIORITY_INIT'][cat['filler']] = 8500
cat['PRIORITY_DONE'] = 0
cat['PRIORITY_DONE'][cat['primary']] = 8000
cat['PRIORITY_DONE'][cat['filler']] = 9000
cat.write('/global/cfs/cdirs/desi/users/rongpu/misc/cosmos_lowz_sample_v0.4_scheme_2.fits')
