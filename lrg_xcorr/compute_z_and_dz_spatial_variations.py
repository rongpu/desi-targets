from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp
from multiprocessing import Pool

nmad = lambda x: 1.4826 * np.median(np.abs(x-np.median(x)))

plot_dir = '/global/cfs/cdirs/desi/users/rongpu/lrg_xcorr/specz_variation/'

nside = 16

min_nobs = 2
# max_ebv = 0.15
# max_stardens = 2500

main_dir = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.1.1/resolve'
tmp = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/original_dr9/dr9_lrg_1.1.1_pzbins_20221204.fits'))
tmp1 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/original_dr9/more/dr9_lrg_1.1.1_pzbins_20221204-weights.fits'))
lrg = hstack([tmp, tmp1], join_type='exact')


# target_bit = 0  # LRG
# fn = '/global/cfs/cdirs/desi/spectro/redux/himalayas/zcatalog/ztile-main-dark-cumulative.fits'
# cat = Table(fitsio.read(fn, columns=['DESI_TARGET']))
# idx = np.where(cat['DESI_TARGET'] & 2**target_bit > 0)[0]
# cat = Table(fitsio.read(fn, rows=idx))
# print(len(cat), len(np.unique(cat['TARGETID'])))

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/himalayas_zcatalog/ztile-main-dark-cumulative-lrg.fits'))
cat.rename_columns(['TARGET_RA', 'TARGET_DEC'], ['RA', 'DEC'])

cat = join(cat, lrg[['TARGETID', 'PIXEL_NOBS_G', 'PIXEL_NOBS_R', 'PIXEL_NOBS_Z', 'pz_bin', 'weight', 'lrg_mask', 'Z_PHOT_MEDIAN']], keys='TARGETID')
cat.rename_column('weight', 'imaging_weight')
print(len(cat))

cat['EFFTIME_LRG'] = 12.15 * cat['TSNR2_LRG']

# Remove FIBERSTATUS!=0 fibers
mask = cat['COADD_FIBERSTATUS']==0
print('FIBERSTATUS   ', np.sum(~mask), np.sum(mask), np.sum(~mask)/len(mask))
cat = cat[mask]

# Remove "no data" fibers
mask = cat['ZWARN'] & 2**9==0
print('No data   ', np.sum(~mask), np.sum(mask), np.sum(~mask)/len(mask))
cat = cat[mask]

# Apply LRG mask
mask = cat['lrg_mask']==0
print('LRG mask', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
cat = cat[mask]

# Require a minimum depth for the cat coadd
min_depth = 800.
mask = cat['EFFTIME_LRG']>min_depth
print('Min depth   ', np.sum(~mask), np.sum(mask), np.sum(~mask)/len(mask))
cat = cat[mask]

# Remove duplicated objects
print(len(cat), len(np.unique(cat['TARGETID'])))
cat.sort('EFFTIME_LRG', reverse=True)
_, idx_keep = np.unique(cat['TARGETID'], return_index=True)
cat = cat[idx_keep]
print(len(cat), len(np.unique(cat['TARGETID'])))

bad_fibers = np.loadtxt('/global/u2/r/rongpu/notebooks/lrg_xcorr/data/bad_fibers_himalayas_20230106.txt', dtype=int)
print(len(bad_fibers))
# bad fibers identified from redshift distributions (https://github.com/desihub/desispec/issues/1946)
additional_bad_fibers = [466, 1008, 1098, 1219, 1251, 2675, 2676, 2677, 2678, 2679, 2680, 3994, 3995, 4349, 4720, 2250, 2251, 2252, 2253, 3038]
bad_fibers = np.unique(np.concatenate([bad_fibers, additional_bad_fibers]))
print(len(bad_fibers))
print(bad_fibers)
mask_bad = np.in1d(cat['FIBER'], bad_fibers)
print('Bad fibers', np.sum(~mask_bad), np.sum(mask_bad), np.sum(mask_bad)/len(mask_bad))
cat = cat[~mask_bad]

mask = (cat['PIXEL_NOBS_G']>=min_nobs) & (cat['PIXEL_NOBS_R']>=min_nobs) & (cat['PIXEL_NOBS_Z']>=min_nobs)
print('NOBS', np.sum(mask), np.sum(~mask), np.sum(mask)/len(mask))
cat = cat[mask]

# # Martin's EBV cut
# mask = cat['EBV']<max_ebv
# print('EBV', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
# cat = cat[mask]

# # Martin's STARDENS cut
# stardens = np.load('/global/cfs/cdirs/desi/users/rongpu/useful/healpix_maps/pixweight-dr7.1-0.22.0_stardens_64_ring.npy')
# stardens_nside = 64
# mask = stardens>=max_stardens
# bad_hp_idx = np.arange(len(stardens))[mask]
# cat_hp_idx = hp.pixelfunc.ang2pix(stardens_nside, cat['RA'], cat['DEC'], lonlat=True, nest=False)
# mask_bad = np.in1d(cat_hp_idx, bad_hp_idx)
# print('STARDENS', np.sum(~mask_bad), np.sum(mask_bad), np.sum(mask_bad)/len(mask_bad))
# cat = cat[~mask_bad]

# # Remove objects classified as STARs
# mask = (cat['SPECTYPE']!='STAR') & (cat['Z']>=0.0003)
# print('Remove objects classified as STARs:', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
# cat = cat[mask]

print(len(cat))

# Redshift quality cut
cat['q'] = cat['ZWARN']==0
cat['q'] &= cat['Z']<1.45
cat['q'] &= cat['DELTACHI2']>15
print(np.sum(~cat['q'])/len(cat))
mask_quality = cat['q'].copy()

print('Quality', np.sum(mask_quality), np.sum(~mask_quality), np.sum(~mask_quality)/len(mask_quality))
cat = cat[mask_quality]

mask_star = (cat['SPECTYPE']=='STAR') | (cat['Z']<0.0003)
print('Stars', np.sum(mask_star), np.sum(~mask_star), np.sum(mask_star)/len(mask_star))


cat_all = cat.copy()

n_processes = 64
npix = hp.nside2npix(nside)


def get_z_variation(pix_idx):

    pix_list = pix_unique[pix_idx]

    hp_table = Table()
    hp_table['HPXPIXEL'] = pix_list
    hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_list, nest=False, lonlat=True)
    bin_str = 'bin_'+str(pz_bin)+'_'
    for col in [bin_str+'z_mean', bin_str+'z_median', bin_str+'z_l68', bin_str+'z_u68', bin_str+'z_l95', bin_str+'z_u95', bin_str+'z_nmad']:
        hp_table[col] = np.zeros(len(hp_table), dtype=float)
    hp_table[bin_str+'n_objects'] = np.zeros(len(hp_table), dtype=int)
    for index in np.arange(len(pix_idx)):
        idx = pixorder[pixcnts[pix_idx[index]]:pixcnts[pix_idx[index]+1]]
        hp_table[bin_str+'z_mean'][index] = np.mean(cat['Z'][idx])
        hp_table[bin_str+'z_median'][index], hp_table[bin_str+'z_l68'][index], hp_table[bin_str+'z_u68'][index], hp_table[bin_str+'z_l95'][index], hp_table[bin_str+'z_u95'][index] =\
            np.percentile(cat['Z'][idx], [50, 16., 84., 2.5, 97.5])
        hp_table[bin_str+'z_nmad'][index] = nmad(cat['Z'][idx])
        hp_table[bin_str+'n_objects'][index] = len(idx)

    return hp_table


def get_dz_variation(pix_idx):

    pix_list = pix_unique[pix_idx]

    hp_table = Table()
    hp_table['HPXPIXEL'] = pix_list
    hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_list, nest=False, lonlat=True)
    bin_str = 'bin_'+str(pz_bin)+'_'
    for col in [bin_str+'dz_mean', bin_str+'dz_median', bin_str+'dz_l68', bin_str+'dz_u68', bin_str+'dz_l95', bin_str+'dz_u95', bin_str+'dz_nmad']:
        hp_table[col] = np.zeros(len(hp_table), dtype=float)
    hp_table[bin_str+'n_objects'] = np.zeros(len(hp_table), dtype=int)
    for index in np.arange(len(pix_idx)):
        idx = pixorder[pixcnts[pix_idx[index]]:pixcnts[pix_idx[index]+1]]
        hp_table[bin_str+'dz_mean'][index] = np.mean(cat['Z'][idx]-cat['Z_PHOT_MEDIAN'][idx])
        hp_table[bin_str+'dz_median'][index], hp_table[bin_str+'dz_l68'][index], hp_table[bin_str+'dz_u68'][index], hp_table[bin_str+'dz_l95'][index], hp_table[bin_str+'dz_u95'][index] =\
            np.percentile(cat['Z'][idx]-cat['Z_PHOT_MEDIAN'][idx], [50, 16., 84., 2.5, 97.5])
        hp_table[bin_str+'dz_nmad'][index] = nmad(cat['Z'][idx]-cat['Z_PHOT_MEDIAN'][idx])
        hp_table[bin_str+'n_objects'][index] = len(idx)

    return hp_table


# z variation
for index, pz_bin in enumerate(range(1, 5)):
    print(pz_bin)
    mask = (cat_all['pz_bin']==pz_bin)
    mask &= ~mask_star
    cat = cat_all[mask].copy()

    pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'], cat['DEC'], lonlat=True)
    pix_unique, pixcnts = np.unique(pix_allobj, return_counts=True)

    pixcnts = np.insert(pixcnts, 0, 0)
    pixcnts = np.cumsum(pixcnts)

    pixorder = np.argsort(pix_allobj)

    pix_idx_split = np.array_split(np.arange(len(pix_unique)), n_processes)

    # start multiple worker processes
    with Pool(processes=n_processes) as pool:
        res = pool.map(get_z_variation, pix_idx_split)

    if index==0:
        hp_table = vstack(res)
    else:
        tmp = vstack(res)
        tmp.remove_columns(['RA', 'DEC'])
        hp_table = join(hp_table, tmp, keys='HPXPIXEL', join_type='outer')

hp_table = hp_table.filled(0)
hp_table.sort('HPXPIXEL')
hp_table.write('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/misc/lrg_pzbins_specz_stats_{}.fits'.format(nside))


# dz variation
for index, pz_bin in enumerate(range(1, 5)):
    print(pz_bin)
    mask = (cat_all['pz_bin']==pz_bin)
    mask &= ~mask_star
    cat = cat_all[mask].copy()

    pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'], cat['DEC'], lonlat=True)
    pix_unique, pixcnts = np.unique(pix_allobj, return_counts=True)

    pixcnts = np.insert(pixcnts, 0, 0)
    pixcnts = np.cumsum(pixcnts)

    pixorder = np.argsort(pix_allobj)

    pix_idx_split = np.array_split(np.arange(len(pix_unique)), n_processes)

    # start multiple worker processes
    with Pool(processes=n_processes) as pool:
        res = pool.map(get_dz_variation, pix_idx_split)

    if index==0:
        hp_table = vstack(res)
    else:
        tmp = vstack(res)
        tmp.remove_columns(['RA', 'DEC'])
        hp_table = join(hp_table, tmp, keys='HPXPIXEL', join_type='outer')

hp_table = hp_table.filled(0)
hp_table.sort('HPXPIXEL')
hp_table.write('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/misc/lrg_pzbins_dz_stats_{}.fits'.format(nside))
