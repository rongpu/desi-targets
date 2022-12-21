from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


wiroc = Table(fitsio.read('/global/cscratch1/sd/dstn/wiro-C/tractor/cus/tractor-custom-034600m04950.fits'))
wiroc.rename_columns(wiroc.colnames, [ii.lower() for ii in wiroc.colnames])
print(len(wiroc))

mask = wiroc['type']!='DUP'
print('DUP', np.sum(mask)/len(mask), np.sum(mask), np.sum(~mask))
wiroc = wiroc[mask]
mask = (wiroc['flux_nb_c']>0) & (wiroc['flux_ivar_nb_c']>0)
print('Positive flux and flux_ivar in WIRO-C', np.sum(mask)/len(mask), np.sum(mask), np.sum(~mask))
wiroc = wiroc[mask]
print(len(wiroc))

wirod = Table(fitsio.read('/global/cscratch1/sd/dstn/wiro-D/tractor/cus/tractor-custom-034600m04950.fits'))
wirod.rename_columns(wirod.colnames, [ii.lower() for ii in wirod.colnames])
print(len(wirod))

mask = wirod['type']!='DUP'
print('DUP', np.sum(mask)/len(mask), np.sum(mask), np.sum(~mask))
wirod = wirod[mask]
mask = (wirod['flux_nb_d']>0) & (wirod['flux_ivar_nb_d']>0)
print('Positive flux and flux_ivar in WIRO-D', np.sum(mask)/len(mask), np.sum(mask), np.sum(~mask))
wirod = wirod[mask]
print(len(wirod))

ramin, ramax, decmin, decmax = np.concatenate([wiroc['ra'], wirod['ra']]).min(), np.concatenate([wiroc['ra'], wirod['ra']]).max(), np.concatenate([wiroc['dec'], wirod['dec']]).min(), np.concatenate([wiroc['dec'], wirod['dec']]).max()
ramin, ramax, decmin, decmax = ramin-0.01, ramax+0.01, decmin-0.01, decmax+0.01

sweep_fns = ['/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/sweep-030m005-040p000.fits', '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/sweep-030m010-040m005.fits']
sweep_stack = []
for sweep_fn in sweep_fns:
    sweep = Table(fitsio.read(sweep_fn, columns=['RA', 'DEC']))
    mask = (sweep['RA']>ramin) & (sweep['RA']<ramax) & (sweep['DEC']>decmin) & (sweep['DEC']<decmax)
    idx = np.where(mask)[0]
    sweep = Table(fitsio.read(sweep_fn, rows=idx))
    print(len(sweep))
    sweep_stack.append(sweep)
cat = vstack(sweep_stack)
print(len(cat))
cat.rename_columns(cat.colnames, [ii.lower() for ii in cat.colnames])

wiroc_colnames_old = ['release', 'brickid', 'brickname', 'objid', 'brick_primary', 'maskbits', 'fitbits', 'type', 'ra', 'dec', 'ra_ivar', 'dec_ivar', 'bx', 'by', 'dchisq', 'flux_nb_c', 'flux_ivar_nb_c', 'fiberflux_nb_c', 'fibertotflux_nb_c', 'apflux_nb_c', 'apflux_resid_nb_c', 'apflux_blobresid_nb_c', 'apflux_ivar_nb_c', 'apflux_masked_nb_c', 'nobs_nb_c', 'rchisq_nb_c', 'fracflux_nb_c', 'fracmasked_nb_c', 'fracin_nb_c', 'ngood_nb_c', 'anymask_nb_c', 'allmask_nb_c', 'psfsize_nb_c', 'psfdepth_nb_c', 'galdepth_nb_c', 'nea_nb_c', 'blob_nea_nb_c', 'sersic', 'sersic_ivar', 'shape_r', 'shape_r_ivar', 'shape_e1', 'shape_e1_ivar', 'shape_e2', 'shape_e2_ivar']
wiroc_colnames_new = ['wiroc_release', 'wiroc_brickid', 'wiroc_brickname', 'wiroc_objid', 'wiroc_brick_primary', 'wiroc_maskbits', 'wiroc_fitbits', 'wiroc_type', 'wiroc_ra', 'wiroc_dec', 'wiroc_ra_ivar', 'wiroc_dec_ivar', 'wiroc_bx', 'wiroc_by', 'wiroc_dchisq', 'wiroc_flux_nb_c', 'wiroc_flux_ivar_nb_c', 'wiroc_fiberflux_nb_c', 'wiroc_fibertotflux_nb_c', 'wiroc_apflux_nb_c', 'wiroc_apflux_resid_nb_c', 'wiroc_apflux_blobresid_nb_c', 'wiroc_apflux_ivar_nb_c', 'wiroc_apflux_masked_nb_c', 'wiroc_nobs_nb_c', 'wiroc_rchisq_nb_c', 'wiroc_fracflux_nb_c', 'wiroc_fracmasked_nb_c', 'wiroc_fracin_nb_c', 'wiroc_ngood_nb_c', 'wiroc_anymask_nb_c', 'wiroc_allmask_nb_c', 'wiroc_psfsize_nb_c', 'wiroc_psfdepth_nb_c', 'wiroc_galdepth_nb_c', 'wiroc_nea_nb_c', 'wiroc_blob_nea_nb_c', 'wiroc_sersic', 'wiroc_sersic_ivar', 'wiroc_shape_r', 'wiroc_shape_r_ivar', 'wiroc_shape_e1', 'wiroc_shape_e1_ivar', 'wiroc_shape_e2', 'wiroc_shape_e2_ivar']
wiroc = wiroc[wiroc_colnames_old]
wiroc.rename_columns(wiroc_colnames_old, wiroc_colnames_new)

wirod_colnames_old = ['release', 'brickid', 'brickname', 'objid', 'brick_primary', 'maskbits', 'fitbits', 'type', 'ra', 'dec', 'ra_ivar', 'dec_ivar', 'bx', 'by', 'dchisq', 'flux_nb_d', 'flux_ivar_nb_d', 'fiberflux_nb_d', 'fibertotflux_nb_d', 'apflux_nb_d', 'apflux_resid_nb_d', 'apflux_blobresid_nb_d', 'apflux_ivar_nb_d', 'apflux_masked_nb_d', 'nobs_nb_d', 'rchisq_nb_d', 'fracflux_nb_d', 'fracmasked_nb_d', 'fracin_nb_d', 'ngood_nb_d', 'anymask_nb_d', 'allmask_nb_d', 'psfsize_nb_d', 'psfdepth_nb_d', 'galdepth_nb_d', 'nea_nb_d', 'blob_nea_nb_d', 'sersic', 'sersic_ivar', 'shape_r', 'shape_r_ivar', 'shape_e1', 'shape_e1_ivar', 'shape_e2', 'shape_e2_ivar']
wirod_colnames_new = ['wirod_release', 'wirod_brickid', 'wirod_brickname', 'wirod_objid', 'wirod_brick_primary', 'wirod_maskbits', 'wirod_fitbits', 'wirod_type', 'wirod_ra', 'wirod_dec', 'wirod_ra_ivar', 'wirod_dec_ivar', 'wirod_bx', 'wirod_by', 'wirod_dchisq', 'wirod_flux_nb_d', 'wirod_flux_ivar_nb_d', 'wirod_fiberflux_nb_d', 'wirod_fibertotflux_nb_d', 'wirod_apflux_nb_d', 'wirod_apflux_resid_nb_d', 'wirod_apflux_blobresid_nb_d', 'wirod_apflux_ivar_nb_d', 'wirod_apflux_masked_nb_d', 'wirod_nobs_nb_d', 'wirod_rchisq_nb_d', 'wirod_fracflux_nb_d', 'wirod_fracmasked_nb_d', 'wirod_fracin_nb_d', 'wirod_ngood_nb_d', 'wirod_anymask_nb_d', 'wirod_allmask_nb_d', 'wirod_psfsize_nb_d', 'wirod_psfdepth_nb_d', 'wirod_galdepth_nb_d', 'wirod_nea_nb_d', 'wirod_blob_nea_nb_d', 'wirod_sersic', 'wirod_sersic_ivar', 'wirod_shape_r', 'wirod_shape_r_ivar', 'wirod_shape_e1', 'wirod_shape_e1_ivar', 'wirod_shape_e2', 'wirod_shape_e2_ivar']
wirod = wirod[wirod_colnames_old]
wirod.rename_columns(wirod_colnames_old, wirod_colnames_new)

idx1, idx2, d2d, d_ra, d_dec = match_coord(wiroc['wiroc_ra'], wiroc['wiroc_dec'], cat['ra'], cat['dec'], search_radius=2., plot_q=True)
cat['match_id'] = -1
cat['match_id'][idx2] = np.arange(len(idx1))
wiroc['match_id'] = -2
wiroc['match_id'][idx1] = np.arange(len(idx1))
cat = join(cat, wiroc, keys='match_id', join_type='left')

idx1, idx2, d2d, d_ra, d_dec = match_coord(wirod['wirod_ra'], wirod['wirod_dec'], cat['ra'], cat['dec'], search_radius=2., plot_q=True)
cat['match_id'] = -1
cat['match_id'][idx2] = np.arange(len(idx1))
wirod['match_id'] = -2
wirod['match_id'][idx1] = np.arange(len(idx1))
cat = join(cat, wirod, keys='match_id', join_type='left')

mask = (~cat['wiroc_ra'].mask) | (~cat['wirod_ra'].mask)
cat = cat[mask]

for col in cat.colnames:
    if not hasattr(cat[col], 'filled'):
        continue
    dtype = type(cat[col][0])
    if dtype==np.str_:
        cat[col] = cat[col].filled('')
    else:
        cat[col] = cat[col].filled(0)

cat.write('/global/cfs/cdirs/desi/users/rongpu/xmm_lae/wiro_cd_dr9_merged.fits')
