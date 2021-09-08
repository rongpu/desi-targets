# Predict the DECam z and GAIA G magnitudes using Tycho-2 and 2MASS photometry

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
import match_coord


tycho2_path = '/global/cfs/cdirs/desi/users/rongpu/useful/Tycho-2.fits'
twomass_path = '/global/cfs/cdirs/desi/users/rongpu/useful/2mass_psc/2mass_psc_j_12.fits'

output_path = '/global/cfs/cdirs/desi/users/rongpu/useful/Tycho-2-reference.fits'

tycho2 = Table(fitsio.read(tycho2_path))
tycho2['idx'] = np.arange(len(tycho2))

twomass = Table(fitsio.read(twomass_path))

mask = np.isfinite(tycho2['RAmdeg']) & np.isfinite(tycho2['DEmdeg'])
mask &= np.isfinite(tycho2['VTmag'])
print(np.sum(~mask)/len(mask))
tycho2 = tycho2[mask]

idx1, idx2, d2d, d_ra, d_dec = match_coord.match_coord(twomass['RAJ2000'], twomass['DEJ2000'], tycho2['RAmdeg'], tycho2['DEmdeg'], priority2=-tycho2['VTmag'], search_radius=5., plot_q=False)
print(len(idx1)/len(twomass))
print(len(idx1)/len(tycho2))

tycho2['Jmag'] = np.nan
tycho2['zguess'] = np.nan
tycho2['ggguess'] = np.nan

twomass = twomass[idx1]
tycho2['Jmag'][idx2] = twomass['Jmag']

coeffs_z = [-0.01835938, -0.68084937, 0.49222576]
coeffs_gg = [0.00445346, -0.07819228, -0.07145574, 0.00278177]
xmin, xmax = -1, 8

x = tycho2['VTmag'][idx2]-twomass['Jmag']

pz = np.poly1d(coeffs_z)
tycho2['zguess'][idx2] = pz(np.clip(x, xmin, xmax)) + tycho2['VTmag'][idx2]

pgg = np.poly1d(coeffs_gg)
tycho2['ggguess'][idx2] = pgg(np.clip(x, xmin, xmax)) + tycho2['VTmag'][idx2]

tycho2.write(output_path)
