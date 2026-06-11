from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord

# south_fns = glob.glob('/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/*.fits')
# north_fns = glob.glob('/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/north/sweep/9.0/*.fits')
# south_fns = [os.path.basename(tmp) for tmp in south_fns]
# north_fns = [os.path.basename(tmp) for tmp in north_fns]
# print(len(south_fns), len(north_fns))

# both_fns = np.intersect1d(south_fns, north_fns)
# print(len(both_fns))

# mask = ['p030-' in tmp for tmp in both_fns]
# both_fns = both_fns[mask]
# print(len(both_fns))

# fns_to_drop = ['sweep-010p030-020p035.fits', 'sweep-330p030-340p035.fits', 'sweep-340p030-350p035.fits', 'sweep-350p030-360p035.fits']
# mask = ~np.in1d(both_fns, fns_to_drop)
# both_fns = both_fns[mask]
# print(len(both_fns))
# print(both_fns)

sweep_fns = ['sweep-100p030-110p035.fits',
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
 'sweep-270p030-280p035.fits']

min_nobs = 1

for sweep_fn in sweep_fns:

    print(sweep_fn)

    sweep_path_north = os.path.join('/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/north/sweep/9.0/', sweep_fn)
    sweep_path_south = os.path.join('/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/', sweep_fn)

    sweep_north = Table(fitsio.read(sweep_path_north))
    sweep_south = Table(fitsio.read(sweep_path_south))

    mask = sweep_north['TYPE']!='DUP'
    mask &= (sweep_north['NOBS_G']>=min_nobs) & (sweep_north['NOBS_R']>=min_nobs) & (sweep_north['NOBS_Z']>=min_nobs)
    print(np.sum(mask), np.sum(mask)/len(mask))
    sweep_north = sweep_north[mask]

    mask = sweep_south['TYPE']!='DUP'
    mask &= (sweep_south['NOBS_G']>=min_nobs) & (sweep_south['NOBS_R']>=min_nobs) & (sweep_south['NOBS_Z']>=min_nobs)
    print(np.sum(mask), np.sum(mask)/len(mask))
    sweep_south = sweep_south[mask]

    # Only keep objects in the overlapping rectangle
    ##################################
    # DOES NOT WORK NEAR THE POLES OR RA=0!!!!
    #################################
    ra1min, ra1max, dec1min, dec1max = sweep_north['RA'].min(), sweep_north['RA'].max(), sweep_north['DEC'].min(), sweep_north['DEC'].max()
    ra1min, ra1max, dec1min, dec1max = ra1min - 10/3600, ra1max + 10/3600, dec1min - 5/3600, dec1max + 5/3600
    ra2min, ra2max, dec2min, dec2max = sweep_south['RA'].min(), sweep_south['RA'].max(), sweep_south['DEC'].min(), sweep_south['DEC'].max()
    ra2min, ra2max, dec2min, dec2max = ra2min - 10/3600, ra1max + 10/3600, dec1min - 5/3600, dec1max + 5/3600

    mask = (sweep_north['RA']>ra2min) & (sweep_north['RA']<ra2max) & (sweep_north['DEC']>dec2min) & (sweep_north['DEC']<dec2max)
    sweep_north = sweep_north[mask]
    mask = (sweep_south['RA']>ra1min) & (sweep_south['RA']<ra1max) & (sweep_south['DEC']>dec1min) & (sweep_south['DEC']<dec1max)
    sweep_south = sweep_south[mask]

    if len(sweep_north)==0 or len(sweep_south)==0:
        continue

    idx1, idx2, d2d, d_ra, d_dec = match_coord(sweep_north['RA'], sweep_north['DEC'], sweep_south['RA'], sweep_south['DEC'], search_radius=.5, plot_q=False)

    sweep_north = sweep_north[idx1]
    sweep_south = sweep_south[idx2]

    sweep_north.write('/global/cscratch1/sd/rongpu/dr9_tests/targets/north_south_overlap/'+sweep_fn.replace('.fits', '-north.fits'))
    sweep_south.write('/global/cscratch1/sd/rongpu/dr9_tests/targets/north_south_overlap/'+sweep_fn.replace('.fits', '-south.fits'))
