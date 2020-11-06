# Select *subsample* of broadly-defined LRGs with photo-z's in DR9 south

# Require NOBS>=2 for grz
# Using DR8 photo-z's
# No extinction correction

from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
import gc
from multiprocessing import Pool

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord

field = 'south'
# field = 'north'

downsample_factor = 16
n_processess = 32

sweep_dir = '/global/cscratch1/sd/adamyers/dr9m-sep26-2020/{}/sweep'.format(field)
sweep_dr8_dir = '/global/project/projectdirs/cosmo/data/legacysurvey/dr8/{}/sweep/8.0'.format(field)
pz_dr8_dir = '/global/project/projectdirs/cosmo/data/legacysurvey/dr8/{}/sweep/8.0-photo-z'.format(field)

output_path = '/Users/rongpu/Documents/Data/desi_lrg_selection/dr9m-sep26-2020/lrg_extended_20201016_ds_{}.fits'.format(field)

sweep_all_path = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))
sweep_fn_list = [os.path.basename(sweep_all_path[ii]) for ii in range(len(sweep_all_path))]

columns = ['TYPE', 'RA', 'DEC', 'EBV', 
          'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2', 
          'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2', 
          'NOBS_G', 'NOBS_R', 'NOBS_Z', 'NOBS_W1', 'NOBS_W2', 
          'FRACFLUX_G', 'FRACFLUX_R', 'FRACFLUX_Z', 'FRACFLUX_W1', 'FRACFLUX_W2', 
          'FRACMASKED_G', 'FRACMASKED_R', 'FRACMASKED_Z', 'FRACIN_G', 'FRACIN_R', 'FRACIN_Z', 
          'WISEMASK_W1', 'WISEMASK_W2', 'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z', 
          'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z', 'PSFDEPTH_W1', 'PSFDEPTH_W2', 
          'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 
          'GAIA_PHOT_G_MEAN_MAG', 'MASKBITS']

pz_columns = ['z_phot_mean', 'z_phot_median', 'z_phot_std', 'z_phot_l68', 'z_phot_u68', 'z_phot_l95', 'z_phot_u95',
              'z_spec', 'survey', 'training']

def trim_sweep(sweep_index):
    
    sweep_fn = sweep_fn_list[sweep_index]
    print(sweep_fn)

    cat = Table(fitsio.read(os.path.join(sweep_dir, sweep_fn),
                columns=['FLUX_Z', 'FIBERFLUX_Z', 'FLUX_IVAR_Z', 'MW_TRANSMISSION_Z']))
    # Apply the (magnitude cut zmag<22.0) OR (fibermag cut zfibermag<22.5)
    mask = (cat['FLUX_Z']>0) & (cat['FLUX_IVAR_Z']>0)
    mask &= (cat['FLUX_Z']/cat['MW_TRANSMISSION_Z'] > 10**(0.4*(22.5-22.0))) | (cat['FIBERFLUX_Z']/cat['MW_TRANSMISSION_Z'] > 10**(0.4*(22.5-22.5)))
    if np.sum(mask)==0:
        return None
    idx = np.where(mask)[0]

    if len(idx)>=downsample_factor:
        idx = np.random.choice(idx, len(idx)//downsample_factor, replace=False)
        idx.sort()
    else:
        return None

    cat = Table(fitsio.read(os.path.join(sweep_dir, sweep_fn), columns=columns, rows=idx))
    cat['TYPE'] = cat['TYPE'].astype(str)

    cat['FLUX_R_EC'] = cat['FLUX_R']/cat['MW_TRANSMISSION_R']
    cat['FLUX_Z_EC'] = cat['FLUX_Z']/cat['MW_TRANSMISSION_Z']
    cat['FLUX_W1_EC'] = cat['FLUX_W1']/cat['MW_TRANSMISSION_W1']

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Quality cuts
        mask = (cat['NOBS_G']>=2) & (cat['NOBS_R']>=2) & (cat['NOBS_Z']>=2)
        # GAIA duplicates (they all have empty photometry values):
        mask &= (cat['TYPE']!='DUP') & (cat['TYPE']!='DUP ')
        
        # # Quality in r: SNR_R > 0 && RFLUX > 0
        # mask &= (cat['FLUX_R']>0) & (cat['FLUX_IVAR_R']>0)

        # Quality in z: SNR_Z > 0 && ZFLUX > 0
        mask &= (cat['FLUX_Z']>0) & (cat['FLUX_IVAR_Z']>0)

        # Quality in W1: FLUX_IVAR_W1 > 0 && W1FLUX > 0
        mask &= (cat['FLUX_W1']>0) & (cat['FLUX_IVAR_W1']>0)

        # None-stellar color: (z-w1) > 0.8*(r-z) - 1.1 => -0.8*r + 1.8*z - W1 > -1.1
        mask_stellar = (cat['FLUX_R_EC']**(-0.8) * cat['FLUX_Z_EC']**1.8 / cat['FLUX_W1_EC'] < 10**(-0.4*(-1.1)))
        # Include non-point sources
        mask_stellar |= ((cat['TYPE']!='PSF') & (cat['TYPE']!='PSF '))

        mask &= mask_stellar

    if np.sum(mask)>0:
        cat = cat[mask]
    else:
        return None

    # Remove unwanted columns
    cat.remove_columns(['FLUX_R_EC', 'FLUX_Z_EC', 'FLUX_W1_EC'])

    # Add DR8 photo-z's
    cat_dr8_path = os.path.join(sweep_dr8_dir, sweep_fn)
    pz_dr8_path = os.path.join(pz_dr8_dir, sweep_fn[:-5]+'-pz.fits')
    
    if os.path.isfile(cat_dr8_path):
        
        cat_dr8 = Table(fitsio.read(cat_dr8_path, columns=['RA', 'DEC', 'TYPE']))
        cat_dr8['TYPE'] = cat_dr8['TYPE'].astype(str)
        pz_dr8 = Table(fitsio.read(pz_dr8_path, columns=pz_columns))
        pz_dr8['survey'] = pz_dr8['survey'].astype(str)
        
        mask = (cat_dr8['TYPE']!='DUP') & (cat_dr8['TYPE']!='DUP ')
        cat_dr8 = cat_dr8[mask]
        pz_dr8 = pz_dr8[mask]

        idx1, idx2, d2d, d_ra, d_dec = match_coord(cat_dr8['RA'], cat_dr8['DEC'], cat['RA'], cat['DEC'], search_radius=1., plot_q=False)

        # Create a padded photo-z table
        tmp = pz_dr8[np.zeros(len(cat), dtype=int)]
        for col in pz_columns:
            if col=='survey':
                tmp[col] = ''
            elif col=='training':
                tmp[col] = False
            else:
                tmp[col] = 0
        tmp[idx2] = pz_dr8[idx1]
        cat = hstack([cat, tmp])

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

    # "Fix" bug in vstack that creates excessively long strings
    cat_stack['TYPE'] = cat_stack['TYPE'].astype('a4')

    print('Final combined catalog:', len(cat_stack))

    cat_stack.write(output_path)

    print(time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
