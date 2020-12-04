from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


north_sweeps = glob.glob('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/north/sweep/9.0/*.fits')
north_sweep_fns = [os.path.basename(tmp) for tmp in north_sweeps]

south_sweeps = glob.glob('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/south/sweep/9.0/*.fits')
south_sweep_fns = [os.path.basename(tmp) for tmp in south_sweeps]

sweep_fns = list(np.intersect1d(north_sweep_fns, south_sweep_fns))
print(sweep_fns)


sweep_fns = ['sweep-010p030-020p035.fits',
'sweep-100p030-110p035.fits',
'sweep-110p030-120p035.fits',
'sweep-120p030-130p035.fits',
'sweep-130p030-140p035.fits',
'sweep-140p030-150p035.fits',
'sweep-150p030-160p035.fits',
'sweep-160p030-170p035.fits',
'sweep-170p030-180p035.fits',
'sweep-180p030-190p035.fits',
'sweep-190p030-200p035.fits',
'sweep-200p030-210p035.fits',
'sweep-210p030-220p035.fits',
'sweep-220p030-230p035.fits',
'sweep-230p030-240p035.fits',
'sweep-240p030-250p035.fits',
'sweep-250p030-260p035.fits',
'sweep-260p030-270p035.fits',
'sweep-270p030-280p035.fits',
'sweep-330p030-340p035.fits',
'sweep-340p030-350p035.fits',
'sweep-350p030-360p035.fits']

sweep_north_all = []
sweep_south_all = []

for sweep_fn in sweep_fns:
    print(sweep_fn)
    sweep_path_north = os.path.join('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/north/sweep/9.0/', sweep_fn)
    sweep_path_south = os.path.join('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/south/sweep/9.0/', sweep_fn)
    
    sweep_north = Table(fitsio.read(sweep_path_north))
    sweep_south = Table(fitsio.read(sweep_path_south))
    
    mask = sweep_north['TYPE']!='DUP'
    mask &= (sweep_north['NOBS_G']>=3) & (sweep_north['NOBS_R']>=3) & (sweep_north['NOBS_Z']>=3)
    print(np.sum(mask), np.sum(mask)/len(mask))
    sweep_north = sweep_north[mask]

    mask = sweep_south['TYPE']!='DUP'
    mask &= (sweep_south['NOBS_G']>=3) & (sweep_south['NOBS_R']>=3) & (sweep_south['NOBS_Z']>=3)
    print(np.sum(mask), np.sum(mask)/len(mask))
    sweep_south = sweep_south[mask]
    
    if len(sweep_north)==0 or len(sweep_south)==0:
        continue
    
    idx1, idx2, d2d, d_ra, d_dec = match_coord(sweep_north['RA'], sweep_north['DEC'], sweep_south['RA'], sweep_south['DEC'], search_radius=.5, plot_q=False)
    
    sweep_north = sweep_north[idx1]
    sweep_south = sweep_south[idx2]
    
    sweep_north_all.append(sweep_north)
    sweep_south_all.append(sweep_south)
    
sweep_north_all = vstack(sweep_north_all)
sweep_south_all = vstack(sweep_south_all)
print(len(sweep_north_all))

sweep_north_all.write('/global/cscratch1/sd/rongpu/dr9_tests/targets/north_vs_south/north.fits')
sweep_south_all.write('/global/cscratch1/sd/rongpu/dr9_tests/targets/north_vs_south/south.fits')
