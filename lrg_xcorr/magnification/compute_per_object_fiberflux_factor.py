from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

from scipy.interpolate import interp1d
from scipy.interpolate import griddata


columns = ['MASKBITS', 'RA', 'DEC', 'TYPE', 'EBV', 'FLUX_Z', 'FIBERFLUX_Z', 'SHAPE_R', 'SHAPE_E1', 'SHAPE_E2', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'SERSIC']

for field in ['north', 'south']:
    
    cat = Table(fitsio.read('/Users/rongpu/Documents/Data/lrg_xcorr/magnification/lrg_magnification_{}.fits'.format(field), columns=columns))

    # axis ratio
    e = np.array(np.sqrt(cat['SHAPE_E1']**2+cat['SHAPE_E2']**2))
    cat['q'] = (1+e)/(1-e)

    cat['ff_ratio'] = -99.  # FIBERFLUX / FLUX
    cat['ff_factor'] = -99.  # Fiberflux multiplicative factor for magnification

    mask_psf = cat['TYPE']=='PSF'
    mask_rex = cat['TYPE']=='REX'
    mask_exp = cat['TYPE']=='EXP'
    mask_dev = cat['TYPE']=='DEV'
    mask_ser = cat['TYPE']=='SER'

    ################################### PSF ###################################

    print('PSF')
    print(np.sum(mask_psf), np.sum(mask_psf)/len(mask_psf))
    cat['ff_factor'][mask_psf] = 1.

    ################################### REX ###################################

    print('REX')
    print(np.sum(mask_rex), np.sum(mask_rex)/len(mask_rex))

    data_rex = np.load('/Users/rongpu/git/desi-targets/lrg_xcorr/magnification/fiberflux_magnification/data/rex.npz')
    f_ratio_interp = interp1d(data_rex['shape_r'], data_rex['ratio'], bounds_error=False, fill_value='extrapolate', kind='quadratic')
    f_ratio = f_ratio_interp(cat['SHAPE_R'][mask_rex])
    f_ratio = np.clip(f_ratio, 0, 1)
    f_factor_interp = interp1d(data_rex['shape_r'], data_rex['f_factor'], bounds_error=False, fill_value='extrapolate', kind='quadratic')
    f_factor = f_factor_interp(cat['SHAPE_R'][mask_rex])
    f_factor = np.clip(f_factor, 0, 1)
    print(np.median(f_factor))

    cat['ff_ratio'][mask_rex] = f_ratio
    cat['ff_factor'][mask_rex] = f_factor

    ################################### EXP ###################################

    print('EXP')
    print(np.sum(mask_exp), np.sum(mask_exp)/len(mask_exp))

    data_exp_ratio = np.load('/Users/rongpu/git/desi-targets/lrg_xcorr/magnification/fiberflux_magnification/data/exp_fiber_ratio.npz')
    def f_ratio_interp_exp(r, q):
        q = np.clip(q, 1.2, 9.95)
        r = np.clip(r, 0., 5.99)
        f = griddata(np.array([data_exp_ratio['shape_r'], data_exp_ratio['q']]).T, data_exp_ratio['ratio'], (r, q), method='cubic')
        f = np.clip(f, 0., 1.)
        return f
    f_ratio = f_ratio_interp_exp(cat['SHAPE_R'][mask_exp], cat['q'][mask_exp])
    mask = np.isnan(f_ratio)
    print('nan', np.sum(mask), np.sum(mask)/np.sum(mask_exp))
    # f_ratio[mask] = np.median(f_ratio[~mask])
    print(np.nanmedian(f_ratio))
    cat['ff_ratio'][mask_exp] = f_ratio

    data_exp_factor = np.load('/Users/rongpu/git/desi-targets/lrg_xcorr/magnification/fiberflux_magnification/data/exp_fiber_factor.npz')
    def f_factor_interp_exp(r, q):
        q = np.clip(q, 1.2, 9.95)
        r = np.clip(r, 0., 5.99)
        f = griddata(np.array([data_exp_factor['shape_r'], data_exp_factor['q']]).T, data_exp_factor['f_factor'], (r, q), method='cubic')
        f = np.clip(f, 0., 1.)
        return f
    f_factor = f_factor_interp_exp(cat['SHAPE_R'][mask_exp], cat['q'][mask_exp])
    mask = np.isnan(f_factor)
    print('nan', np.sum(mask), np.sum(mask)/np.sum(mask_exp))
    f_factor[mask] = np.median(f_factor[~mask])
    print(np.median(f_factor))
    cat['ff_factor'][mask_exp] = f_factor

    ################################### DEV ###################################

    print('DEV')
    print(np.sum(mask_dev), np.sum(mask_dev)/len(mask_dev))

    data_dev_ratio = np.load('/Users/rongpu/git/desi-targets/lrg_xcorr/magnification/fiberflux_magnification/data/dev_fiber_ratio.npz')
    def f_ratio_interp_dev(r, q):
        q = np.clip(q, 1.03, 9.98)
        r = np.clip(r, 0., 9.99)
        f = griddata(np.array([data_dev_ratio['shape_r'], data_dev_ratio['q']]).T, data_dev_ratio['ratio'], (r, q), method='cubic')
        f = np.clip(f, 0., 1.)
        return f
    f_ratio = f_ratio_interp_dev(cat['SHAPE_R'][mask_dev], cat['q'][mask_dev])
    mask = np.isnan(f_ratio)
    print('nan', np.sum(mask), np.sum(mask)/np.sum(mask_dev))
    # f_ratio[mask] = np.median(f_ratio[~mask])
    print(np.nanmedian(f_ratio))
    cat['ff_ratio'][mask_dev] = f_ratio

    data_dev_factor = np.load('/Users/rongpu/git/desi-targets/lrg_xcorr/magnification/fiberflux_magnification/data/dev_fiber_factor.npz')
    def f_factor_interp_dev(r, q):
        q = np.clip(q, 1.03, 9.98)
        r = np.clip(r, 0., 9.99)
        f = griddata(np.array([data_dev_factor['shape_r'], data_dev_factor['q']]).T, data_dev_factor['f_factor'], (r, q), method='cubic')
        f = np.clip(f, 0., 1.)
        return f
    f_factor = f_factor_interp_dev(cat['SHAPE_R'][mask_dev], cat['q'][mask_dev])
    mask = np.isnan(f_factor)
    print('nan', np.sum(mask), np.sum(mask)/np.sum(mask_dev))
    f_factor[mask] = np.median(f_factor[~mask])
    print(np.median(f_factor))
    cat['ff_factor'][mask_dev] = f_factor

    ################################### SER ###################################

    print('SER')
    print(np.sum(mask_ser), np.sum(mask_ser)/len(mask_ser))

    print('SER-EXP')
    mask_ser_exp = mask_ser & (cat['SERSIC']<2.5)

    f_ratio = f_ratio_interp_exp(cat['SHAPE_R'][mask_ser_exp], cat['q'][mask_ser_exp])    
    mask = np.isnan(f_ratio)
    print('nan', np.sum(mask), np.sum(mask)/np.sum(mask_ser))
    # f_ratio[mask] = np.median(f_ratio[~mask])
    print(np.nanmedian(f_ratio))
    cat['ff_ratio'][mask_ser_exp] = f_ratio

    f_factor = f_factor_interp_exp(cat['SHAPE_R'][mask_ser_exp], cat['q'][mask_ser_exp])    
    mask = np.isnan(f_factor)
    print('nan', np.sum(mask), np.sum(mask)/np.sum(mask_ser))
    f_factor[mask] = np.median(f_factor[~mask])
    print(np.median(f_factor))
    cat['ff_factor'][mask_ser_exp] = f_factor

    print('SER-DEV')
    mask_ser_dev = mask_ser & (cat['SERSIC']>=2.5)

    f_ratio = f_ratio_interp_dev(cat['SHAPE_R'][mask_ser_dev], cat['q'][mask_ser_dev])    
    mask = np.isnan(f_ratio)
    print('nan', np.sum(mask), np.sum(mask)/np.sum(mask_ser))
    # f_ratio[mask] = np.median(f_ratio[~mask])
    print(np.nanmedian(f_ratio))
    cat['ff_ratio'][mask_ser_dev] = f_ratio

    f_factor = f_factor_interp_dev(cat['SHAPE_R'][mask_ser_dev], cat['q'][mask_ser_dev])    
    mask = np.isnan(f_factor)
    print('nan', np.sum(mask), np.sum(mask)/np.sum(mask_ser))
    f_factor[mask] = np.median(f_factor[~mask])
    print(np.median(f_factor))
    cat['ff_factor'][mask_ser_dev] = f_factor

    if np.sum(cat['ff_factor']<0):
        raise ValueError

    cat.remove_columns(columns)
    cat.write('/Users/rongpu/Documents/Data/lrg_xcorr/magnification/lrg_magnification_{}_fiberflux.fits'.format(field), overwrite=True)
