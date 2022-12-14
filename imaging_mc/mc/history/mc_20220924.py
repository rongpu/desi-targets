from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp
from multiprocessing import Pool
from scipy.interpolate import RectBivariateSpline

sys.path.append(os.path.expanduser('~/git/desi-targets/dr9/create_target_catalogs/main/'))
from select_desi_targets import select_elg_simplified


n_processes = 128
repeats = 256
n_randoms_catalogs = 32

# from https://www.legacysurvey.org/dr9/catalogs/#galactic-extinction-coefficients
ext_coeffs = {'u': 3.995, 'g': 3.214, 'r': 2.165, 'i': 1.592, 'z': 1.211, 'y': 1.064}

# from psfex_fwhm_vs_nea.ipynb
fwhm_scaling = {'g': 0.963, 'r': 0.960, 'z': 0.952}


get_nea = {}
for band in ['g', 'r', 'z']:
    nea_path = '/global/cfs/cdirs/desi/users/rongpu/imaging_mc/nea/nea_vs_fwhm_{}_1024.fits'.format(band)
    nea = Table(fitsio.read(nea_path))
    nea_arr = np.array(nea['nea']).T
    hdr = fitsio.read_header(nea_path, ext=1)
    shape_r_grid = np.arange(hdr['R_MIN'], hdr['R_MAX']+hdr['R_DELTA'], hdr['R_DELTA'])
    fwhm_grid = np.array(nea['fwhm_bin'])
    get_nea[band] = RectBivariateSpline(shape_r_grid, fwhm_grid, nea_arr).ev


def quicksim(truth, cat):
    '''
    truth columns: flux_grz, fiberflux_g, shape_r
    cat columns: psfsize_grz, psfdepth_grz, ebv
    '''
    nea = {}
    for band in ['g', 'r', 'z']:
        nea[band] = np.array(get_nea[band](truth['shape_r'], fwhm_scaling[band]*cat['psfsize_'+band]), dtype='float32')
    # nea = Table(nea)

    flux_err = {}
    for band in ['g', 'r', 'z']:
        pix_ivar = cat['psfdepth_'+band] * 4 * np.pi * (cat['psfsize_'+band]/2.3548)**2  # pixel-level inverse variance per arcsec^2
        flux_err[band] = np.sqrt(nea[band]/pix_ivar)

    sim = Table()
    for band in ['g', 'r', 'z']:
        sim['flux_'+band] = np.array(truth['flux_{}_ec'.format(band)] / 10**(0.4*ext_coeffs[band]*cat['ebv']) + np.random.randn(len(cat)) * flux_err[band], dtype='float32')
    for band in ['g']:
        fiberflux_ratio = truth['fiberflux_{}_ec'.format(band)] / truth['flux_{}_ec'.format(band)]
        sim['fiberflux_'+band] = np.array(truth['fiberflux_{}_ec'.format(band)] / 10**(0.4*ext_coeffs[band]*cat['ebv']) + np.random.randn(len(cat)) * fiberflux_ratio*flux_err[band], dtype='float32')

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for band in ['g', 'r', 'z']:
            sim[band+'mag'] = 22.5 - 2.5*np.log10(sim['flux_'+band]) - ext_coeffs[band] * cat['ebv']
        for band in ['g']:
            sim[band+'fibermag'] = 22.5 - 2.5*np.log10(sim['fiberflux_'+band]) - ext_coeffs[band] * cat['ebv']

    sim = sim[['gmag', 'rmag', 'zmag', 'gfibermag']]

    return sim


def elgsim(foo):

    if len(cat)<=len(truth):
        replace = False
    else:
        replace = True

    np.random.seed((os.getpid() * int(time.time()*1000)) % 123456789)

    idx = np.random.choice(len(truth), size=len(cat), replace=replace)
    sim = quicksim(truth[idx], cat)

    elglop, elgvlo = select_elg_simplified(sim)

    return elglop, elgvlo


truth = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/imaging_mc/subsets/cosmos_truth_clean.fits'))
print('truth', len(truth))

# extinction-corrected fluxes
for band in ['g', 'r', 'z']:
    truth['flux_{}_ec'.format(band)] = truth['flux_'+band]*10**(0.4*ext_coeffs[band]*truth['ebv'])
for band in ['g']:
    truth['fiberflux_{}_ec'.format(band)] = truth['fiberflux_'+band]*10**(0.4*ext_coeffs[band]*truth['ebv'])

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    truth['gfibermag'] = 22.5 - 2.5*np.log10(truth['fiberflux_g']) - ext_coeffs['g'] * truth['ebv']

mask = truth['gfibermag']<25.1
truth = truth[mask]
print('truth', len(truth))

# Only keep necessary columns to save memory
truth = truth[['flux_g_ec', 'flux_r_ec', 'flux_z_ec', 'fiberflux_g_ec', 'shape_r']]

# cat_north = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/resolve/combined/pixmap_north_nside_128_minobs_1_maskbits__elgmask_v1.fits'))
# cat_south = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/resolve/combined/pixmap_south_nside_128_minobs_1_maskbits__elgmask_v1.fits'))
# cat_north['PHOTSYS'] = 'N'
# cat_south['PHOTSYS'] = 'S'
# mask = (cat_north['DEC']>32.375)
# cat_north = cat_north[mask]
# mask = ~np.in1d(cat_south['HPXPIXEL'], cat_north['HPXPIXEL'])
# cat = vstack([cat_north, cat_south[mask]])
# print('cat', len(cat))

columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'EBV', 'PHOTSYS']

randoms_paths = sorted(glob.glob('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-[0-9]*.fits'))
randoms_paths = randoms_paths[:n_randoms_catalogs]

for randoms_path in randoms_paths:
    print(randoms_path)

    output_path = '/global/cfs/cdirs/desi/users/rongpu/imaging_mc/mc/mc_elg_'+os.path.basename(randoms_path).replace('.fits', '-elgmask_v1.fits')
    # output_path = '/pscratch/sd/r/rongpu/imaging_mc/mc_elg_'+os.path.basename(randoms_path).replace('.fits', '-elgmask_v1.fits')
    if os.path.isfile(output_path):
        continue

    cat = Table(fitsio.read(randoms_path, columns=columns))
    cat1 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/desi_mask/randoms/elgmask_v1/'+os.path.basename(randoms_path).replace('.fits', '-elgmask_v1.fits')))
    cat = hstack([cat, cat1])
    min_nobs = 1
    mask = (cat['NOBS_G']>=min_nobs) & (cat['NOBS_R']>=min_nobs) & (cat['NOBS_Z']>=min_nobs)
    mask &= (cat['PSFDEPTH_G']>0) & (cat['PSFDEPTH_R']>0) & (cat['PSFDEPTH_Z']>0)
    mask &= (cat['PSFSIZE_G']>0) & (cat['PSFSIZE_R']>0) & (cat['PSFSIZE_Z']>0)
    cat = cat[mask]
    print('cat', len(cat))
    mask = cat['elg_mask']==0
    cat = cat[mask]
    print('cat', len(cat))

    cat.rename_columns(cat.colnames, [ii.lower() for ii in cat.colnames])

    elglop_count, elgvlo_count = np.zeros(len(cat), dtype=int), np.zeros(len(cat), dtype=int)
    counter = repeats
    while counter>0:
        with Pool(processes=n_processes) as pool:
            res = pool.map(elgsim, np.zeros(np.minimum(counter, n_processes)))
        counter = np.maximum(counter-n_processes, 0)
        res = np.array(res, dtype=int)
        elglop_count += np.sum(res[:, 0], axis=0)
        elgvlo_count += np.sum(res[:, 1], axis=0)
        del res

    mc = Table()
    mc['ra'] = cat['ra']
    mc['dec'] = cat['dec']
    mc['photsys'] = cat['photsys']
    mc['elglop'] = elglop_count
    mc['elgvlo'] = elgvlo_count
    mc.write(output_path, overwrite=True)

################################################# Healpix density map ###########################################################

def healpix_stats(pix_idx):

    pix_list = pix_unique[pix_idx]

    hp_table = Table()
    hp_table['HPXPIXEL'] = pix_list
    hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_list, nest=False, lonlat=True)
    hp_table['elglop'] = np.zeros(len(hp_table), dtype=int)
    hp_table['elgvlo'] = np.zeros(len(hp_table), dtype=int)

    for index in np.arange(len(pix_idx)):

        idx = pixorder[pixcnts[pix_idx[index]]:pixcnts[pix_idx[index]+1]]
        hp_table['elglop'][index] = np.sum(mc['elglop'][idx])
        hp_table['elgvlo'][index] = np.sum(mc['elgvlo'][idx])

    return hp_table


nside = 128

mc_stack = []
for randoms_path in randoms_paths:
    output_path = '/global/cfs/cdirs/desi/users/rongpu/imaging_mc/mc/mc_elg_'+os.path.basename(randoms_path).replace('.fits', '-elgmask_v1.fits')
    # output_path = '/pscratch/sd/r/rongpu/imaging_mc/mc_elg_'+os.path.basename(randoms_path).replace('.fits', '-elgmask_v1.fits')
    mc_stack.append(Table(fitsio.read(output_path)))
mc = vstack(mc_stack)

pix_allobj = hp.pixelfunc.ang2pix(nside, mc['ra'], mc['dec'], nest=False, lonlat=True)
pix_unique, pix_count = np.unique(pix_allobj, return_counts=True)
pixcnts = pix_count.copy()
pixcnts = np.insert(pixcnts, 0, 0)
pixcnts = np.cumsum(pixcnts)
pixorder = np.argsort(pix_allobj)
# split among the processors
pix_idx_split = np.array_split(np.arange(len(pix_unique)), n_processes)
# start multiple worker processes
with Pool(processes=n_processes) as pool:
    res = pool.map(healpix_stats, pix_idx_split)
hp_table = vstack(res)
hp_table.sort('HPXPIXEL')
hp_table['n_randoms'] = pix_count

hp_table.write('/global/cfs/cdirs/desi/users/rongpu/imaging_mc/mc/mc_elg_randoms_elgmask_v1_healpix_{}.fits'.format(nside), overwrite=True)
# hp_table.write('/pscratch/sd/r/rongpu/imaging_mc/mc_elg_randoms_elgmask_v1_healpix_{}.fits'.format(nside), overwrite=True)

