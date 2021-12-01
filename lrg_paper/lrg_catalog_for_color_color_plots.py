from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
import gc
from astropy.io import fits
from multiprocessing import Pool

field = 'south'
downsample_factor = 16
n_processess = 4

sweep_dir = '/global/project/projectdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0'.format(field)
pz_dir = '/global/project/projectdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0-photo-z'.format(field)
output_path = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/paper/lrg_extended_20211122_{}.fits'.format(field)

# sweep_all_path = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))
# sweep_fn_list = [os.path.basename(sweep_all_path[ii]) for ii in range(len(sweep_all_path))]

sweep_fn_list = ['sweep-130p000-140p005.fits','sweep-140p000-150p005.fits','sweep-150p000-160p005.fits','sweep-160p000-170p005.fits','sweep-170p000-180p005.fits','sweep-180p000-190p005.fits','sweep-190p000-200p005.fits','sweep-200p000-210p005.fits','sweep-210p000-220p005.fits','sweep-220p000-230p005.fits','sweep-230p000-240p005.fits','sweep-240p000-250p005.fits','sweep-130p005-140p010.fits','sweep-140p005-150p010.fits','sweep-150p005-160p010.fits','sweep-160p005-170p010.fits','sweep-170p005-180p010.fits','sweep-180p005-190p010.fits','sweep-190p005-200p010.fits','sweep-200p005-210p010.fits','sweep-210p005-220p010.fits','sweep-220p005-230p010.fits','sweep-230p005-240p010.fits','sweep-240p005-250p010.fits','sweep-130p010-140p015.fits','sweep-140p010-150p015.fits','sweep-150p010-160p015.fits','sweep-160p010-170p015.fits','sweep-170p010-180p015.fits','sweep-180p010-190p015.fits','sweep-190p010-200p015.fits','sweep-200p010-210p015.fits','sweep-210p010-220p015.fits','sweep-220p010-230p015.fits','sweep-230p010-240p015.fits','sweep-240p010-250p015.fits']

columns = ['RELEASE', 'BRICKID', 'OBJID', 'TYPE', 'RA', 'DEC', 'EBV',
           'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2',
           'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2',
           'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2',
           'FIBERFLUX_Z', 'FIBERTOTFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG',
           'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'WISEMASK_W1']

np.random.seed(852)

def trim_sweep(sweep_index):
    
    sweep_fn = sweep_fn_list[sweep_index]
    print(sweep_fn)

    cat = Table(fitsio.read(os.path.join(sweep_dir, sweep_fn),
            columns=['FLUX_Z', 'FIBERFLUX_Z', 'FLUX_IVAR_Z', 'MW_TRANSMISSION_Z']))

    # Apply the (magnitude cut zmag<21.) OR (fibermag cut zfibermag<22)
    zmag = 22.5 - 2.5 * np.log10((cat['FLUX_Z'] / cat['MW_TRANSMISSION_Z']).clip(1e-7))
    zfibermag = 22.5 - 2.5 * np.log10((cat['FIBERFLUX_Z'] / cat['MW_TRANSMISSION_Z']).clip(1e-7))
    mask = (cat['FLUX_Z']>0) & (cat['FLUX_IVAR_Z']>0)
    mask &= (zmag<21.) | (zfibermag<22.)
    if np.sum(mask)==0:
        return None
    idx = np.where(mask)[0]

    if len(idx)>=downsample_factor:
        idx = np.random.choice(idx, len(idx)//downsample_factor, replace=False)
        idx.sort()
    else:
        return None

    cat = Table(fitsio.read(os.path.join(sweep_dir, sweep_fn), columns=columns, rows=idx))

    pz_path = os.path.join(pz_dir, sweep_fn[:-5]+'-pz.fits')
    pz = Table(fitsio.read(pz_path, rows=idx))

    cat = hstack([cat, pz])

    gmag = 22.5 - 2.5 * np.log10((cat['FLUX_G'] / cat['MW_TRANSMISSION_G']).clip(1e-7))
    # ADM safe as these fluxes are set to > 0 in notinLRG_mask.
    rmag = 22.5 - 2.5 * np.log10((cat['FLUX_R'] / cat['MW_TRANSMISSION_R']).clip(1e-7))
    zmag = 22.5 - 2.5 * np.log10((cat['FLUX_Z'] / cat['MW_TRANSMISSION_Z']).clip(1e-7))
    w1mag = 22.5 - 2.5 * np.log10((cat['FLUX_W1'] / cat['MW_TRANSMISSION_W1']).clip(1e-7))
    zfibermag = 22.5 - 2.5 * np.log10((cat['FIBERFLUX_Z'] / cat['MW_TRANSMISSION_Z']).clip(1e-7))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Quality cuts
        mask = (cat['NOBS_G']>=2) & (cat['NOBS_R']>=2) & (cat['NOBS_Z']>=2)
        mask &= (cat['TYPE']!='DUP') & (cat['TYPE']!='DUP ')
        
        # # Quality in r: SNR_R > 0 && RFLUX > 0
        # mask &= (cat['FLUX_R_EC']>0) & (cat['FLUX_IVAR_R']>0)

        # Quality in z: SNR_Z > 0 && ZFLUX > 0
        mask &= (cat['FLUX_Z']>0) & (cat['FLUX_IVAR_Z']>0)

        # Quality in W1: FLUX_IVAR_W1 > 0 && W1FLUX > 0
        mask &= (cat['FLUX_W1']>0) & (cat['FLUX_IVAR_W1']>0)

        # None-stellar color: (z-w1) > 0.8*(r-z) - 1.0
        mask_stellar = zmag - w1mag > 0.8 * (rmag - zmag) - 1.0  # non-stellar cut
        # Include non-point sources
        mask_stellar |= ((cat['TYPE']!='PSF') & (cat['TYPE']!='PSF '))

        mask &= mask_stellar

    if np.sum(mask)>0:
        cat = cat[mask]
    else:
        return None

    # clear cache
    gc.collect()

    return cat


if __name__ == '__main__':

    print('Start!')

    time_start = time.time()

    # start multiple worker processes
    with Pool(processes=n_processess) as pool:
        res = pool.map(trim_sweep, range(len(sweep_fn_list)))

    # Remove None elements from the list
    for index in range(len(res)-1, -1, -1):
        if res[index] is None:
            res.pop(index)

    cat_stack = vstack(res)

    print('Final combined catalog:', len(cat_stack))

    cat_stack.write(output_path)

    print(time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
