from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from multiprocessing import Pool
from scipy.interpolate import RectBivariateSpline

sys.path.append(os.path.expanduser('~/git/desi-targets/dr9/create_target_catalogs/main/'))
from select_desi_targets import select_elg_simplified


n_processes = 128

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
    truth columns: flux_grz, fiberflux_g, shape_r, ebv
    cat columns: psfsize_grz, psfdepth_grz, ebv
    '''

    # from psfex_fwhm_vs_nea.ipynb
    fwhm_scaling = {'g': 0.963, 'r': 0.960, 'z': 0.952}

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
        sim['flux_'+band] = np.array(truth['flux_'+band] + np.random.randn(len(flux_err[band])) * flux_err[band], dtype='float32')
    for band in ['g', 'r', 'z']:
        fiberflux_ratio = truth['fiberflux_'+band] / truth['flux_'+band]
        sim['fiberflux_'+band] = np.array(truth['fiberflux_'+band] + np.random.randn(len(flux_err[band])) * fiberflux_ratio*flux_err[band], dtype='float32')

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        sim['gmag'] = 22.5 - 2.5*np.log10(sim['flux_g']) - 3.214 * cat['ebv']
        sim['rmag'] = 22.5 - 2.5*np.log10(sim['flux_r']) - 2.165 * cat['ebv']
        sim['zmag'] = 22.5 - 2.5*np.log10(sim['flux_z']) - 1.211 * cat['ebv']
        sim['gfibermag'] = 22.5 - 2.5*np.log10(sim['fiberflux_g']) - 3.214 * cat['ebv']
        # sim['rfibermag'] = 22.5 - 2.5*np.log10(sim['fiberflux_r']) - 2.165 * cat['ebv']
        # sim['zfibermag'] = 22.5 - 2.5*np.log10(sim['fiberflux_z']) - 1.211 * cat['ebv']

    sim = sim[['gmag', 'rmag', 'zmag', 'gfibermag']]

    return sim


def elgsim():

    if len(cat)<=len(truth):
        replace = False
    else:
        replace = True
    idx = np.random.choice(len(truth), size=len(cat), replace=replace)
    sim = quicksim(truth[idx], cat)

    elglop, elgvlo = select_elg_simplified(sim)

    return elglop, elgvlo


# def count_elg(repeats=1, random_seed=None):

truth = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/imaging_mc/subsets/cosmos_truth_clean.fits'))
print('truth', len(truth))

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    truth['gfibermag'] = 22.5 - 2.5*np.log10(truth['fiberflux_g']) - 3.214 * truth['ebv']

mask = truth['gfibermag']<25.1
truth = truth[mask]
print('truth', len(truth))

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/resolve/combined/pixmap_south_nside_256_minobs_1_maskbits__elgmask_v1.fits'))
print('cat', len(cat))
cat.rename_columns(cat.colnames, [ii.lower() for ii in cat.colnames])

# np.random.seed(random_seed)
repeats = 1000

with Pool(processes=n_processes) as pool:
    res = pool.starmap(elgsim, [[]]* repeats)
res = np.array(res, dtype=int)
elglop_count, elgvlo_count = np.sum(res[:, 0], axis=0), np.sum(res[:, 1], axis=0)

tt = Table()
tt['elglop'] = elglop_count
tt['elgvlop'] = elgvlo_count
tt.write('/global/u2/r/rongpu/temp/data/mc_elg.fits')

