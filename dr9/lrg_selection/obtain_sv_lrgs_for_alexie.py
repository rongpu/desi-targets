# Create DR9 SV LRG catalog for Alexie

from __future__ import division, print_function
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
# import matplotlib.pyplot as plt
import sys, os, glob, time, warnings, gc
from multiprocessing import Pool

n_processes = 32
field = 'south'
output_path = '/global/cscratch1/sd/rongpu/share/dr9_lrg/dr9_sv_lrg_radec_only_{}.fits'.format(field)

sweep_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/{}/sweep/9.0'.format(field)

columns = ['RA', 'DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1',
        'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FIBERFLUX_Z',
        'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS',
        'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1']

columns_to_keep = ['RA', 'DEC', 'LRG_OPT', 'LRG_IR', 'LRG_SV_OPT', 'LRG_SV_IR']

# for sweep_index in range(len(sweep_fn_list)):
def get_lrgs(sweep_index):

    sweep_fn = sweep_fn_list[sweep_index]
    print(sweep_fn)
    cat = Table(fitsio.read(os.path.join(sweep_dir, sweep_fn), columns=columns))

    # Apply masks
    maskbits = [1, 12, 13]
    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    print('{:} ({:.1f}%) objects removed by maskbits'.format(np.sum(~mask_clean), np.sum(~mask_clean)/len(mask_clean)*100))
    cat = cat[mask_clean]
    print(len(cat))

    # Basic quality cuts
    mask = (cat['FLUX_R']>0) & (cat['FLUX_IVAR_R']>0)
    mask &= (cat['FLUX_Z']>0) & (cat['FIBERFLUX_Z']>0) & (cat['FLUX_IVAR_Z']>0)
    mask &= (cat['FLUX_W1']>0) & (cat['FLUX_IVAR_W1']>0)
    # print(np.sum(~mask)/len(mask), np.sum(~mask))
    cat = cat[mask]
    print(len(cat))

    # Require >=1 visits in g,r and z
    mask = (cat['NOBS_G']>=1) & (cat['NOBS_R']>=1) & (cat['NOBS_Z']>=1)
    cat = cat[mask]
    print(len(cat))

    if np.sum(mask)==0:
        return None

    # Replace negative fluxes by 1e-7 (40 mag)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_G']/cat['MW_TRANSMISSION_G'], 1e-7, None))
        rmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_R']/cat['MW_TRANSMISSION_R'], 1e-7, None))
        zmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_Z']/cat['MW_TRANSMISSION_Z'], 1e-7, None))
        zmag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_Z']/cat['MW_TRANSMISSION_Z'], 1e-7, None))
        w1mag = 22.5 - 2.5*np.log10(np.clip(cat['FLUX_W1']/cat['MW_TRANSMISSION_W1'], 1e-7, None))
        zfibermag = 22.5 - 2.5*np.log10(np.clip(cat['FIBERFLUX_Z']/cat['MW_TRANSMISSION_Z'], 1e-7, None))

    if field=='south':

        ############ LRG_OPT ############
        lrg_mask = zmag - w1mag > 0.8 * (rmag-zmag) - 0.6  # non-stellar cut
        lrg_mask &= (zfibermag < 21.5)                     # faint limit
        mask_red = (gmag - w1mag > 2.6) & (gmag - rmag > 1.4)  # low-z cut
        mask_red |= (rmag-w1mag) > 1.8                         # ignore low-z cut for faint objects
        lrg_mask &= mask_red
        lrg_mask &= rmag - zmag > (zmag - 16.83) * 0.45       # sliding optical cut
        lrg_mask &= rmag - zmag > (zmag - 13.80) * 0.19       # low-z sliding optical cut
        lrg_opt = lrg_mask.copy()

        ############ LRG_IR ############
        lrg_mask = zmag - w1mag > 0.8 * (rmag-zmag) - 0.6  # non-stellar cut
        lrg_mask &= (zfibermag < 21.5)                     # faint limit
        lrg_mask &= (rmag - w1mag > 1.1)                   # Low-z cut
        lrg_mask &= rmag - w1mag > (w1mag - 17.22) * 1.8   # sliding IR cut
        lrg_mask &= rmag - w1mag > w1mag - 16.37           # low-z sliding IR cut
        lrg_ir = lrg_mask.copy()

        ############ LRG_SV_OPT ############
        lrg_mask = zmag - w1mag > 0.8 * (rmag-zmag) - 0.8  # non-stellar cut
        lrg_mask &= ((zmag < 21.0) | (zfibermag < 22.0))   # faint limit
        mask_red = (gmag - w1mag > 2.5) & (gmag - rmag > 1.3)  # low-z cut
        mask_red |= (rmag-w1mag) > 1.7                         # ignore low-z cut for faint objects
        lrg_mask &= mask_red
        # straight cut for low-z:
        lrg_mask_lowz = zmag < 20.2
        lrg_mask_lowz &= rmag - zmag > (zmag - 17.20) * 0.45
        lrg_mask_lowz &= rmag - zmag > (zmag - 14.17) * 0.19
        # curved sliding cut for high-z:
        lrg_mask_highz = zmag >= 20.2
        lrg_mask_highz &= (((zmag - 23.18) / 1.3)**2 + (rmag - zmag + 2.5)**2 > 4.48**2)
        lrg_mask &= (lrg_mask_lowz | lrg_mask_highz)
        lrg_sv_opt = lrg_mask.copy()

        ############ LRG_SV_IR ############
        lrg_mask = zmag - w1mag > 0.8 * (rmag-zmag) - 0.8  # non-stellar cut
        lrg_mask &= ((zmag < 21.0) | (zfibermag < 22.0))   # faint limit
        lrg_mask &= (rmag - w1mag > 1.0)                   # Low-z cut
        lrg_mask_slide = rmag - w1mag > (w1mag - 17.48) * 1.8  # sliding IR cut
        lrg_mask_slide |= (rmag - w1mag > 3.1)
        lrg_mask &= lrg_mask_slide
        lrg_sv_ir = lrg_mask.copy()

    elif field=='north':

        ############ LRG_OPT ############
        lrg_mask = zmag - w1mag > 0.8 * (rmag-zmag) - 0.6  # non-stellar cut
        lrg_mask &= (zfibermag < 21.5)                     # faint limit
        mask_red = (gmag - w1mag > 2.67) & (gmag - rmag > 1.45)  # low-z cut
        mask_red |= (rmag-w1mag) > 1.85                          # ignore low-z cut for faint objects
        lrg_mask &= mask_red
        lrg_mask &= rmag - zmag > (zmag - 16.79) * 0.45       # sliding optical cut
        lrg_mask &= rmag - zmag > (zmag - 13.76) * 0.19       # low-z sliding optical cut
        lrg_opt = lrg_mask.copy()

        ############ LRG_IR ############
        lrg_mask = zmag - w1mag > 0.8 * (rmag-zmag) - 0.6  # non-stellar cut
        lrg_mask &= (zfibermag < 21.5)                     # faint limit
        lrg_mask &= (rmag - w1mag > 1.13)                  # Low-z cut
        lrg_mask &= rmag - w1mag > (w1mag - 17.18) * 1.8   # sliding IR cut
        lrg_mask &= rmag - w1mag > w1mag - 16.33           # low-z sliding IR cut
        lrg_ir = lrg_mask.copy()

        ############ LRG_SV_OPT ############
        lrg_mask = zmag - w1mag > 0.8 * (rmag-zmag) - 0.8  # non-stellar cut
        lrg_mask &= ((zmag < 21.0) | (zfibermag < 22.0))   # faint limit
        mask_red = (gmag - w1mag > 2.57) & (gmag - rmag > 1.35)  # low-z cut
        mask_red |= (rmag-w1mag) > 1.75                          # ignore low-z cut for faint objects
        lrg_mask &= mask_red
        # straight cut for low-z:
        lrg_mask_lowz = zmag < 20.2
        lrg_mask_lowz &= rmag - zmag > (zmag - 17.17) * 0.45
        lrg_mask_lowz &= rmag - zmag > (zmag - 14.14) * 0.19
        # curved sliding cut for high-z:
        lrg_mask_highz = zmag >= 20.2
        lrg_mask_highz &= (((zmag - 23.15) / 1.3)**2 + (rmag - zmag + 2.5)**2 > 4.48**2)
        lrg_mask &= (lrg_mask_lowz | lrg_mask_highz)
        lrg_sv_opt = lrg_mask.copy()

        ############ LRG_SV_IR ############
        lrg_mask = zmag - w1mag > 0.8 * (rmag-zmag) - 0.8  # non-stellar cut
        lrg_mask &= ((zmag < 21.0) | (zfibermag < 22.0))   # faint limit
        lrg_mask &= (rmag - w1mag > 1.03)                  # Low-z cut
        lrg_mask_slide = rmag - w1mag > (w1mag - 17.44) * 1.8  # sliding IR cut
        lrg_mask_slide |= (rmag - w1mag > 3.1)
        lrg_mask &= lrg_mask_slide
        lrg_sv_ir = lrg_mask.copy()

    lrg_all = lrg_opt | lrg_ir | lrg_sv_opt | lrg_sv_ir
    if np.sum(lrg_all)==0:
        return None

    print('{:.3f}%'.format(np.sum(lrg_all)/len(lrg_all)*100), np.sum(lrg_all))
    cat = cat[lrg_all]

    cat['LRG_OPT'] = np.array(lrg_opt[lrg_all], dtype=np.int16)
    cat['LRG_IR'] = np.array(lrg_ir[lrg_all], dtype=np.int16)
    cat['LRG_SV_OPT'] = np.array(lrg_sv_opt[lrg_all], dtype=np.int16)
    cat['LRG_SV_IR'] = np.array(lrg_sv_ir[lrg_all], dtype=np.int16)

    # Remove unwanted columns
    cat = cat[columns_to_keep]

    # clear cache
    gc.collect()

    return cat


if __name__ == '__main__':

    print('Start!')

    time_start = time.time()

    sweep_all_path = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))
    sweep_fn_list = [os.path.basename(tmp) for tmp in sweep_all_path]

    # start multiple worker processes
    with Pool(processes=n_processes) as pool:
        res = pool.map(get_lrgs, range(len(sweep_fn_list)))

    # Remove None elements from the list
    for index in range(len(res)-1, -1, -1):
        if res[index] is None:
            res.pop(index)

    cat_stack = vstack(res)

    print('Final combined catalog:', len(cat_stack))

    cat_stack.write(output_path, overwrite=True)

    print(time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
