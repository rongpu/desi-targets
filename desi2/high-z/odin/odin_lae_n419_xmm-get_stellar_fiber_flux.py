from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/desi-examples/misc/misc'))
from estimate_fiber_flux_from_stars import get_stellar_flux
cat = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/users/dstn/ODIN/xmm-N419/tractor-xmm-N419-hsc-forced.fits'))
ffcat = get_stellar_flux(cat)

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    ffcat['stellar_gmag'] = 22.5 - 2.5*np.log10(ffcat['star_flux_g'])
    ffcat['stellar_rmag'] = 22.5 - 2.5*np.log10(ffcat['star_flux_r'])
    ffcat['stellar_imag'] = 22.5 - 2.5*np.log10(ffcat['star_flux_i'])
    ffcat['stellar_zmag'] = 22.5 - 2.5*np.log10(ffcat['star_flux_z'])

ffcat.write('/global/cfs/cdirs/desi/users/rongpu/xmm_lae/test/misc/tractor-xmm-N419-hsc-forced_stellar_fiber_flux.fits', overwrite=True)
