# Based on lrg_catalog_for_color_color_plots.py

from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
import gc
from astropy.io import fits
from multiprocessing import Pool

field = 'south'
downsample_factor = 8
n_processess = 18

sweep_dir = '/global/project/projectdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0'.format(field)
pz_dir = '/global/project/projectdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0-photo-z'.format(field)
output_path = '/global/cfs/cdirs/desi/users/rongpu/data/desi2/parent_catalog_20220324.fits'

# sweep_all_path = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))
# sweep_fn_list = [os.path.basename(sweep_all_path[ii]) for ii in range(len(sweep_all_path))]

sweep_fn_list = ['sweep-130p000-140p005.fits', 'sweep-140p000-150p005.fits', 'sweep-150p000-160p005.fits', 'sweep-160p000-170p005.fits', 'sweep-170p000-180p005.fits', 'sweep-180p000-190p005.fits', 'sweep-190p000-200p005.fits', 'sweep-200p000-210p005.fits', 'sweep-210p000-220p005.fits', 'sweep-220p000-230p005.fits', 'sweep-230p000-240p005.fits', 'sweep-240p000-250p005.fits', 'sweep-130p005-140p010.fits', 'sweep-140p005-150p010.fits', 'sweep-150p005-160p010.fits', 'sweep-160p005-170p010.fits', 'sweep-170p005-180p010.fits', 'sweep-180p005-190p010.fits', 'sweep-190p005-200p010.fits', 'sweep-200p005-210p010.fits', 'sweep-210p005-220p010.fits', 'sweep-220p005-230p010.fits', 'sweep-230p005-240p010.fits', 'sweep-240p005-250p010.fits', 'sweep-130p010-140p015.fits', 'sweep-140p010-150p015.fits', 'sweep-150p010-160p015.fits', 'sweep-160p010-170p015.fits', 'sweep-170p010-180p015.fits', 'sweep-180p010-190p015.fits', 'sweep-190p010-200p015.fits', 'sweep-200p010-210p015.fits', 'sweep-210p010-220p015.fits', 'sweep-220p010-230p015.fits', 'sweep-230p010-240p015.fits', 'sweep-240p010-250p015.fits']

columns = ['RELEASE', 'BRICKID', 'OBJID', 'TYPE', 'RA', 'DEC', 'EBV',
           'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2',
           'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2',
           'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 'GAIA_PHOT_G_MEAN_MAG',
           'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'WISEMASK_W1']

np.random.seed(852)


def trim_sweep(sweep_index):

    sweep_fn = sweep_fn_list[sweep_index]
    print(sweep_fn)

    cat = Table(fitsio.read(os.path.join(sweep_dir, sweep_fn),
                columns=['FLUX_Z', 'FIBERFLUX_Z', 'FLUX_IVAR_Z', 'FLUX_R', 'FIBERFLUX_R', 'FLUX_IVAR_R', 'EBV']))

    # Apply the magnitude cuts
    zmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_Z']*10**(0.4*1.211*cat['EBV']), 1e-7, None))
    zfibermag = 22.5 - 2.5*np.log10(np.clip(cat['FIBERFLUX_Z']*10**(0.4*1.211*cat['EBV']), 1e-7, None))
    rmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_R']*10**(0.4*2.165*cat['EBV']), 1e-7, None))
    rfibermag = 22.5 - 2.5*np.log10(np.clip(cat['FIBERFLUX_R']*10**(0.4*2.165*cat['EBV']), 1e-7, None))
    # mask = (cat['FLUX_Z']>0) & (cat['FLUX_IVAR_Z']>0)
    # mask &= (cat['FLUX_R']>0) & (cat['FLUX_IVAR_R']>0)
    mask = (cat['FLUX_IVAR_Z']>0) & (cat['FLUX_IVAR_R']>0)
    mask &= (zmag<21.) | (zfibermag<22.5) | (rmag<22.) | (rfibermag<23.5)
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
    pz.remove_columns(['RELEASE', 'BRICKID', 'OBJID'])

    cat = hstack([cat, pz])

    gmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_G']*10**(0.4*3.214*cat['EBV']), 1e-7, None))
    rmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_R']*10**(0.4*2.165*cat['EBV']), 1e-7, None))
    zmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_Z']*10**(0.4*1.211*cat['EBV']), 1e-7, None))
    w1mag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_W1']*10**(0.4*0.184*cat['EBV']), 1e-7, None))
    w2mag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_W2']*10**(0.4*0.113*cat['EBV']), 1e-7, None))
    rfibermag = 22.5 - 2.5*np.log10(np.clip(cat['FIBERFLUX_R']*10**(0.4*2.165*cat['EBV']), 1e-7, None))
    zfibermag = 22.5 - 2.5*np.log10(np.clip(cat['FIBERFLUX_Z']*10**(0.4*1.211*cat['EBV']), 1e-7, None))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Quality cuts
        mask = (cat['NOBS_G']>=2) & (cat['NOBS_R']>=2) & (cat['NOBS_Z']>=2)
        mask &= (cat['TYPE']!='DUP') & (cat['TYPE']!='DUP ')

        # # Quality in r: SNR_R > 0 && RFLUX > 0
        # mask &= (cat['FLUX_R_EC']>0) & (cat['FLUX_IVAR_R']>0)

        # Quality in z: SNR_Z > 0 && ZFLUX > 0
        mask &= (cat['FLUX_Z']>0) & (cat['FLUX_IVAR_Z']>0)

        # # Quality in W1: FLUX_IVAR_W1 > 0 && W1FLUX > 0
        # mask &= (cat['FLUX_W1']>0) & (cat['FLUX_IVAR_W1']>0)

        # # None-stellar color: (z-w1) > 0.8*(r-z) - 1.0
        # mask_stellar = zmag - w1mag > 0.8 * (rmag - zmag) - 1.1  # non-stellar cut
        # # Include non-point sources
        # mask_stellar |= ((cat['TYPE']!='PSF') & (cat['TYPE']!='PSF '))
        # mask &= mask_stellar

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

    cat_stack.write(output_path, overwrite=True)

    print(time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
