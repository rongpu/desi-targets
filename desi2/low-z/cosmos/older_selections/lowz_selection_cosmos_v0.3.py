from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio


# COSMOS field center
ra, dec = 150.11917, 2.20583
ramin, ramax, decmin, decmax = ra - 2, ra + 2, dec - 2, dec + 2

# sweep_fns = ['sweep-140p000-150p005.fits', 'sweep-150p000-160p005.fits']
# print(len(sweep_fns))

# columns = ['RELEASE', 'BRICKID', 'OBJID', 'RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'TYPE', 'EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG']

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
mask_quality &= (cat['FLUX_R'] > 0)   # ADM quality in r.
mask_quality &= (cat['FLUX_Z'] > 0) & (cat['FIBERFLUX_Z'] > 0)   # ADM quality in z.
mask_quality &= (cat['FLUX_W1'] > 0)  # ADM quality in W1.
mask_quality &= (cat['GAIA_PHOT_G_MEAN_MAG'] == 0) | (cat['GAIA_PHOT_G_MEAN_MAG'] > 18)  # remove bright GAIA sources
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

############################ Selection cuts ############################
gmag = 22.5 - 2.5 * np.log10((cat['FLUX_G']*10**(0.4*3.214*cat['EBV'])).clip(1e-7))
# ADM safe as these fluxes are set to > 0 in notinLRG_mask.
rmag = 22.5 - 2.5 * np.log10((cat['FLUX_R']*10**(0.4*2.165*cat['EBV'])).clip(1e-7))
zmag = 22.5 - 2.5 * np.log10((cat['FLUX_Z']*10**(0.4*1.211*cat['EBV'])).clip(1e-7))
w1mag = 22.5 - 2.5 * np.log10((cat['FLUX_W1']*10**(0.4*0.184*cat['EBV'])).clip(1e-7))
zfibermag = 22.5 - 2.5 * np.log10((cat['FIBERFLUX_Z']*10**(0.4*1.211*cat['EBV'])).clip(1e-7))

mask = zfibermag<21.6
mask &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # LRG non-stellar cut
mask &= (gmag - w1mag > 2.9) | (rmag - w1mag > 1.4)  # low-z cuts
print(np.sum(mask), np.sum(mask)/len(cat))
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

#################### Exclude SV3 LRG and SV3 BGS_BRIGHT+BGS_FAINT targets ######################
sv3_dir = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/0.57.0'

sv3lrg = []
for field in ['north', 'south']:
    tmp = Table(fitsio.read(os.path.join(sv3_dir, 'dr9_sv3_lrg_{}_0.57.0_basic.fits'.format(field)), columns=['TARGETID']))
    sv3lrg.append(hstack([tmp]))
sv3lrg = vstack(sv3lrg)
print(len(sv3lrg))

sv3bgs = []
for field in ['north', 'south']:
    tmp = Table(fitsio.read(os.path.join(sv3_dir, 'dr9_sv3_bgs_any_{}_0.57.0_basic.fits'.format(field)), columns=['TARGETID', 'SV3_BGS_TARGET']))
    sv3bgs.append(hstack([tmp]))
sv3bgs = vstack(sv3bgs)
print(len(sv3bgs))
# Select BGS Bright and BGS Faint targets
mask = (sv3bgs['SV3_BGS_TARGET'] & 2**1 > 0) | (sv3bgs['SV3_BGS_TARGET'] & 2**0 > 0)
sv3bgs = sv3bgs[mask]
print(len(sv3bgs))

cat['sv3lrg'] = np.in1d(cat['TARGETID'], sv3lrg['TARGETID'])
cat['sv3bgs'] = np.in1d(cat['TARGETID'], sv3bgs['TARGETID'])

cat['is_target'] = (~cat['sv3lrg']) & (~cat['sv3bgs'])

cat.write('/global/cfs/cdirs/desi/users/rongpu/misc/cosmos_lowz_sample_v0.3.fits')

