# Select DESI targets based on the "intrinsic" magnitudes and colors are corrected for EBV (based on DESI stellar spectra) and zero-point offsets (based on Gaia-LS comparison)

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, join
import fitsio
import healpy as hp

nside = 128

# Load delta_gr and delta_rz maps; RING orderinrg
maps = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/v0/desi_std/maps/delta_gr_map_all_{}.fits'.format(nside)))
maps.rename_columns(['EBV', 'n_star'], ['EBV_SFD', 'n_star_gr'])
# maps['EBV_GR'] = maps['delta_gr_hlmean']/1.049
grmaps = maps.copy()

maps = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/v0/desi_std/maps/delta_rz_map_all_{}.fits'.format(nside)))
maps.remove_columns(['EBV'])
maps.rename_columns(['n_star'], ['n_star_rz'])
# maps['EBV_RZ'] = maps['delta_rz_hlmean']/0.954
rzmaps = maps.copy()

maps = join(rzmaps, grmaps, keys='HPXPIXEL', join_type='inner')
print(len(maps), len(maps)/len(grmaps), len(maps)/len(rzmaps))

# Load ZP offset maps; RING orderinrg
dr9_south_offsets = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/misc/gaia_xp_dr9_south_offset_maps_{}.fits'.format(nside)))
dr9_north_offsets = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/misc/gaia_xp_dr9_north_offset_maps_{}.fits'.format(nside)))

mask = (dr9_north_offsets['DEC']>32.375) & (dr9_north_offsets['RA']>90) & (dr9_north_offsets['RA']<300)
dr9_north_offsets = dr9_north_offsets[mask]

mask = ~np.in1d(dr9_south_offsets['HPXPIXEL'], dr9_north_offsets['HPXPIXEL'])
dr9_south_offsets = dr9_south_offsets[mask]

dr9_offsets = vstack([dr9_south_offsets, dr9_north_offsets])
dr9_offsets.sort('HPXPIXEL')

maps = join(maps, dr9_offsets, keys='HPXPIXEL', join_type='inner')
print(len(maps))

mask = np.isfinite(maps['gmag_diff_median']) & np.isfinite(maps['rmag_diff_median']) & np.isfinite(maps['zmag_diff_median'])
maps = maps[mask]
print(np.sum(mask)/len(mask), len(maps))

maps['EBV_DESI'] = (maps['delta_gr_hlmean']-(maps['gmag_diff_median']-maps['rmag_diff_median']))/1.049


def select_lrg(cat, field='south'):
    '''
    columns = ['OBJID', 'BRICKID', 'RELEASE', 'RA', 'DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FIBERFLUX_Z', 'FIBERTOTFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG', 'EBV', 'MASKBITS', 'NOBS_G', 'NOBS_R', 'NOBS_Z']
    fn = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/sweep-000p000-010p005.fits'
    cat = Table(fitsio.read(fn, columns=columns))
    mask_lrg = select_lrg(cat)
    '''

    cat = cat.copy()
    cat.rename_columns(cat.colnames, [ii.upper() for ii in cat.colnames])

    mask_quality = np.full(len(cat), True)

    mask_quality &= (cat['FLUX_IVAR_R'] > 0) & (cat['FLUX_R'] > 0)   # ADM quality in r.
    mask_quality &= (cat['FLUX_IVAR_Z'] > 0) & (cat['FLUX_Z'] > 0) & (cat['FIBERFLUX_Z'] > 0)   # ADM quality in z.
    mask_quality &= (cat['FLUX_IVAR_W1'] > 0) & (cat['FLUX_W1'] > 0)  # ADM quality in W1.

    mask_quality &= (cat['GAIA_PHOT_G_MEAN_MAG'] == 0) | (cat['GAIA_PHOT_G_MEAN_MAG'] > 18)  # remove bright GAIA sources

    # ADM remove stars with zfibertot < 17.5 that are missing from GAIA.
    mask_quality &= cat['FIBERTOTFLUX_Z'] < 10**(-0.4*(17.5-22.5))

    # ADM observed in every band.
    mask_quality &= (cat['NOBS_G'] > 0) & (cat['NOBS_R'] > 0) & (cat['NOBS_Z'] > 0)

    # Apply masks
    maskbits = [1, 12, 13]
    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    # print(np.sum(~mask_clean)/len(mask_clean))
    mask_quality &= mask_clean

    cat['HPXPIXEL'] = hp.ang2pix(nside, cat['RA'], cat['DEC'], nest=False, lonlat=True)
    cat = join(cat, maps[['HPXPIXEL', 'EBV_DESI', 'gmag_diff_median', 'rmag_diff_median', 'zmag_diff_median']], keys='HPXPIXEL', join_type='inner')

    # gmag = 22.5 - 2.5 * np.log10((cat['FLUX_G'] / cat['MW_TRANSMISSION_G']).clip(1e-7))
    # # ADM safe as these fluxes are set to > 0 in notinLRG_mask.
    # rmag = 22.5 - 2.5 * np.log10((cat['FLUX_R'] / cat['MW_TRANSMISSION_R']).clip(1e-7))
    # zmag = 22.5 - 2.5 * np.log10((cat['FLUX_Z'] / cat['MW_TRANSMISSION_Z']).clip(1e-7))
    # w1mag = 22.5 - 2.5 * np.log10((cat['FLUX_W1'] / cat['MW_TRANSMISSION_W1']).clip(1e-7))
    # zfibermag = 22.5 - 2.5 * np.log10((cat['FIBERFLUX_Z'] / cat['MW_TRANSMISSION_Z']).clip(1e-7))
    gmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_G']*10**(0.4*3.214*cat['EBV_DESI']), 1e-7, None)) - cat['gmag_diff_median']
    rmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_R']*10**(0.4*2.165*cat['EBV_DESI']), 1e-7, None)) - cat['rmag_diff_median']
    zmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_Z']*10**(0.4*1.211*cat['EBV_DESI']), 1e-7, None)) - cat['zmag_diff_median']
    w1mag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_W1']*10**(0.4*0.184*cat['EBV_DESI']), 1e-7, None))
    zfibermag = 22.5 - 2.5*np.log10(np.clip(cat['FIBERFLUX_Z']*10**(0.4*1.211*cat['EBV_DESI']), 1e-7, None)) - cat['zmag_diff_median']

    mask_lrg = mask_quality.copy()

    if field=='south':
        mask_lrg &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
        mask_lrg &= zfibermag < 21.6                   # faint limit
        mask_lrg &= (gmag - w1mag > 2.9) | (rmag - w1mag > 1.8)  # low-z cuts
        mask_lrg &= (
            ((rmag - w1mag > (w1mag - 17.14) * 1.8)
             & (rmag - w1mag > (w1mag - 16.33) * 1.))
            | (rmag - w1mag > 3.3)
        )  # double sliding cuts and high-z extension
    else:
        mask_lrg &= zmag - w1mag > 0.8 * (rmag - zmag) - 0.6  # non-stellar cut
        mask_lrg &= zfibermag < 21.61                   # faint limit
        mask_lrg &= (gmag - w1mag > 2.97) | (rmag - w1mag > 1.8)  # low-z cuts
        mask_lrg &= (
            ((rmag - w1mag > (w1mag - 17.13) * 1.83)
             & (rmag - w1mag > (w1mag - 16.31) * 1.))
            | (rmag - w1mag > 3.4)
        )  # double sliding cuts and high-z extension

    return mask_lrg


def select_elg(cat):
    '''
    columns = ['OBJID', 'BRICKID', 'RELEASE', 'RA', 'DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FIBERFLUX_G', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'EBV', 'MASKBITS', 'NOBS_G', 'NOBS_R', 'NOBS_Z']
    fn = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/sweep-000p000-010p005.fits'
    cat = Table(fitsio.read(fn, columns=columns))
    mask_elglop, mask_elgvlo = select_elg(cat)
    '''

    cat = cat.copy()
    cat.rename_columns(cat.colnames, [ii.upper() for ii in cat.colnames])

    mask_quality = np.full(len(cat), True)

    mask_quality &= (cat['FLUX_IVAR_G'] > 0) & (cat['FLUX_G'] > 0) & (cat['FIBERFLUX_G'] > 0)
    mask_quality &= (cat['FLUX_IVAR_R'] > 0) & (cat['FLUX_R'] > 0)
    mask_quality &= (cat['FLUX_IVAR_Z'] > 0) & (cat['FLUX_Z'] > 0)

    # ADM observed in every band.
    mask_quality &= (cat['NOBS_G'] > 0) & (cat['NOBS_R'] > 0) & (cat['NOBS_Z'] > 0)

    # Apply masks
    maskbits = [1, 12, 13]
    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    # print(np.sum(~mask_clean)/len(mask_clean))
    mask_quality &= mask_clean

    cat['HPXPIXEL'] = hp.ang2pix(nside, cat['RA'], cat['DEC'], nest=False, lonlat=True)
    cat = join(cat, maps[['HPXPIXEL', 'EBV_DESI', 'gmag_diff_median', 'rmag_diff_median', 'zmag_diff_median']], keys='HPXPIXEL', join_type='inner')

    # gmag = 22.5 - 2.5 * np.log10((cat['FLUX_G'] / cat['MW_TRANSMISSION_G']).clip(1e-7))
    # # ADM safe as these fluxes are set to > 0 in notinLRG_mask.
    # rmag = 22.5 - 2.5 * np.log10((cat['FLUX_R'] / cat['MW_TRANSMISSION_R']).clip(1e-7))
    # zmag = 22.5 - 2.5 * np.log10((cat['FLUX_Z'] / cat['MW_TRANSMISSION_Z']).clip(1e-7))
    # w1mag = 22.5 - 2.5 * np.log10((cat['FLUX_W1'] / cat['MW_TRANSMISSION_W1']).clip(1e-7))
    # zfibermag = 22.5 - 2.5 * np.log10((cat['FIBERFLUX_Z'] / cat['MW_TRANSMISSION_Z']).clip(1e-7))
    gmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_G']*10**(0.4*3.214*cat['EBV_DESI']), 1e-7, None)) - cat['gmag_diff_median']
    rmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_R']*10**(0.4*2.165*cat['EBV_DESI']), 1e-7, None)) - cat['rmag_diff_median']
    zmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_Z']*10**(0.4*1.211*cat['EBV_DESI']), 1e-7, None)) - cat['zmag_diff_median']
    gfibermag = 22.5 - 2.5 * np.log10(np.clip(cat['FIBERFLUX_G']*10**(0.4*3.214*cat['EBV_DESI']), 1e-7, None)) - cat['gmag_diff_median']

    mask_elglop = mask_quality.copy()

    mask_elglop &= gmag > 20                       # bright cut.
    mask_elglop &= rmag - zmag > 0.15                  # blue cut.
    mask_elglop &= gfibermag < 24.1  # faint cut.
    mask_elglop &= gmag - rmag < 0.5*(rmag - zmag) + 0.1  # remove stars, low-z galaxies.

    mask_elgvlo = mask_elglop.copy()

    # ADM low-priority OII flux cut.
    mask_elgvlo &= gmag - rmag < -1.2*(rmag - zmag) + 1.6
    mask_elgvlo &= gmag - rmag >= -1.2*(rmag - zmag) + 1.3

    # ADM high-priority OII flux cut.
    mask_elglop &= gmag - rmag < -1.2*(rmag - zmag) + 1.3

    return mask_elglop, mask_elgvlo
