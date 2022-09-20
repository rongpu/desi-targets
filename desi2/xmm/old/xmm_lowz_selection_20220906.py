# This is the old version that was checked in but not observed
# Create target catalog for the XMM low-z secondary program

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

parent_path = '/global/cfs/cdirs/desi/users/rongpu/data/desi2/xmm_lowz/targets/xmm_lowz_targets.fits'
prognum = 2
fa_path = '/global/cfs/cdirs/desi/users/rongpu/data/desi2/xmm_lowz/targets/tertiary-targets-{}.fits'.format(str(prognum).zfill(4))

# XMM field center
ra, dec = 36.448, -4.601
ramin, ramax, decmin, decmax = ra - 2, ra + 2, dec - 2, dec + 2

rosette_ra, rosette_dec = 36.45, -4.60  # for fiber assignment and "in_xmm"

# ############################ Start with sweep catalogs ############################

# sweep_fns = ['sweep-030m005-040p000.fits', 'sweep-030m010-040m005.fits']
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

# cat.write('/global/cscratch1/sd/rongpu/temp/sweep_xmm_secondary.fits', overwrite=True)

############################ Basic quality cuts ############################

cat = Table.read('/global/cscratch1/sd/rongpu/temp/sweep_xmm_secondary.fits')
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
sky1 = SkyCoord(rosette_ra*u.degree, rosette_dec*u.degree, frame='icrs')
sky2 = SkyCoord(cat['RA']*u.degree, cat['DEC']*u.degree, frame='icrs')
mask = np.array(sky2.separation(sky1).to(u.degree))<stats_search_radius
cat['in_xmm'] = mask.copy()
print(np.sum(cat['in_xmm']))

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

principal = cat['zfibermag']<21.6
principal &= (cat['TYPE']!='PSF') | (cat['zmag'] - cat['w1mag'] > 0.8 * (cat['rmag'] - cat['zmag']) - 1.1)  # stellar rejection
print(np.sum(principal), np.sum(principal)/len(cat))

filler = (cat['zfibermag']>=21.6) & (cat['zfibermag']<22.4)
print(np.sum(filler), np.sum(filler)/len(cat))
mask_highz = ((cat['gmag']-cat['rmag']) < 1.2 * (cat['rmag']-cat['zmag'])) | (cat['rmag']-cat['zmag']>1.3)  # high-z cuts
filler_hip = filler & mask_highz
filler_lop = filler & (~mask_highz)
print(np.sum(filler), np.sum(filler_hip)/len(cat))
print(np.sum(filler), np.sum(filler_lop)/len(cat))

cat['principal'] = principal.copy()
cat['filler_hip'] = filler_hip.copy()
cat['filler_lop'] = filler_lop.copy()

mask = cat['principal'] | cat['filler_hip'] | cat['filler_lop']
cat = cat[mask]
print(len(cat))

################### Flag existing DESI redshifts ######################

t1 = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv1-dark.fits'))
t1['bright'] = False
mask = (t1['TARGET_RA']>ramin) & (t1['TARGET_RA']<ramax) & (t1['TARGET_DEC']>decmin) & (t1['TARGET_DEC']<decmax)
t1 = t1[mask]
t2 = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv3-dark.fits'))
t2['bright'] = False
mask = (t2['TARGET_RA']>ramin) & (t2['TARGET_RA']<ramax) & (t2['TARGET_DEC']>decmin) & (t2['TARGET_DEC']<decmax)
t2 = t2[mask]

t3 = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv1-bright.fits'))
t3['bright'] = True
mask = (t3['TARGET_RA']>ramin) & (t3['TARGET_RA']<ramax) & (t3['TARGET_DEC']>decmin) & (t3['TARGET_DEC']<decmax)
t3 = t3[mask]
t4 = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/fuji/zcatalog/zpix-sv3-bright.fits'))
t4['bright'] = True
mask = (t4['TARGET_RA']>ramin) & (t4['TARGET_RA']<ramax) & (t4['TARGET_DEC']>decmin) & (t4['TARGET_DEC']<decmax)
t4 = t4[mask]

# t5 = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/guadalupe/zcatalog/zpix-main-bright.fits'))
# mask = (t5['TARGET_RA']>ramin) & (t5['TARGET_RA']<ramax) & (t5['TARGET_DEC']>decmin) & (t5['TARGET_DEC']<decmax)
# t5 = t5[mask]
# t5['bright'] = False

# obs = vstack([t1, t2, t3, t3, t4, t5], join_type='inner')
obs = vstack([t1, t2, t3, t4], join_type='inner')
print(len(obs), len(np.unique(obs['TARGETID'])))

obs['EFFTIME_LRG'] = 12.15 * obs['TSNR2_LRG']

# Remove duplicates keeping the highest EFFTIME_LRG
obs.sort('EFFTIME_LRG', reverse=True)
_, idx_keep = np.unique(obs['TARGETID'], return_index=True)
obs = obs[idx_keep]
print(len(obs), len(np.unique(obs['TARGETID'])))

obs.rename_columns(['TARGET_RA', 'TARGET_DEC'], ['RA', 'DEC'])
obs.write('/global/cscratch1/sd/rongpu/temp/fuji_in_xmm.fits', overwrite=True)

obs = Table(fitsio.read('/global/cscratch1/sd/rongpu/temp/fuji_in_xmm.fits'))

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

print(np.sum(cat['is_target'] & cat['in_xmm']))

# cat.remove_columns(['in_xmm', 'gmag', 'rmag', 'zmag', 'w1mag', 'w2mag', 'rfibermag', 'zfibermag'])

cat['TARGETID_TERTIARY'] = np.full(len(cat), 0, dtype=np.int64)
mask = cat['is_target'].copy()
cat['TARGETID_TERTIARY'][mask] = encode_targetid(release=8888, brickid=prognum, objid=np.arange(np.sum(mask)))
cat.write(parent_path)

print(len(cat))
print()

###################### Produce the final catalog for fiber assignment ######################

cat = Table(fitsio.read(parent_path))

mask = cat['is_target'].copy()
cat = cat[mask]
print(len(cat))

cat['PMRA'] = 0.
cat['PMDEC'] = 0.
cat['REF_EPOCH'] = 2015.5
init_str = ' ' * np.max([len('PRINCIPAL'), len('FILLER_HIP'), len('FILLER_LOP')])
cat['TERTIARY_TARGET'] = init_str
cat['TERTIARY_TARGET'][cat['principal']] = 'PRINCIPAL'
cat['TERTIARY_TARGET'][cat['filler_hip']] = 'FILLER_HIP'
cat['TERTIARY_TARGET'][cat['filler_lop']] = 'FILLER_LOP'
if np.sum(cat['TERTIARY_TARGET']==init_str)>0:
    raise ValueError

cat['CHECKER'] = 'RZ'

# cat = cat[['TARGETID_TERTIARY', 'RA', 'DEC', 'PMRA', 'PMDEC', 'REF_EPOCH', 'TERTIARY_TARGET', 'CHECKER']]
cat.rename_column('TARGETID', 'TARGETID_ORIGINAL')
cat.rename_column('TARGETID_TERTIARY', 'TARGETID')

# Sanity check
tmp = encode_targetid(release=8888, brickid=prognum, objid=np.arange(len(cat)))
if not np.all(tmp==cat['TARGETID']):
    raise ValueError

# Add required header keywords
hdul = fitsio.FITS(fa_path, mode='rw', clobber=True)
for ext in range(2):
    hdr = fitsio.read_header(parent_path, ext=ext)
    if ext==0:
        hdul.write(None, header=hdr)
    elif ext==1:
        data = np.array(cat.copy())
        hdr['FAPRGRM'] = 'xmmlowz'
        hdr['OBSCONDS'] = 'DARK'
        hdr['SBPROF'] = 'ELG'
        hdr['GOALTIME'] = 1200
        hdul.write(data, header=hdr, extname='TARGETS')
hdul.close()

