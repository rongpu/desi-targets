from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/xmm_lae/wiro_cd_dr9_merged.fits'))
print(len(cat))

# Stellar fiber flux estimates
ffcat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/xmm_lae/misc/wiro_cd_dr9_merged_stellar_fiber_flux.fits'))
cat['stellar_gmag'] = ffcat['stellar_gmag']

mask = cat['stellar_gmag']>20
cat = cat[mask]
print(len(cat))

# Not correcting for extinction
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    cat['gmag'] = 22.5 - 2.5*np.log10(cat['flux_g'])
    cat['rmag'] = 22.5 - 2.5*np.log10(cat['flux_r'])
    cat['zmag'] = 22.5 - 2.5*np.log10(cat['flux_z'])
    cat['nb_cmag'] = 22.5 - 2.5*np.log10(cat['wiroc_flux_nb_c'])
    cat['nb_dmag'] = 22.5 - 2.5*np.log10(cat['wirod_flux_nb_d'])

cmag_synth = cat['gmag'] - 0.15 + (2.25/2)*(cat['gmag']-cat['rmag'])
dmag_synth = cat['gmag'] + 0.10 + (1.20/2)*(cat['gmag']-cat['rmag'])

# mask = cat['gmag']-cat['rmag']<1.2
# plt.figure(figsize=(8, 8))
# plt.plot(np.arange(100), np.arange(100), color='C3')
# plt.plot(cmag_synth[mask], cat['nb_cmag'][mask], '.', ms=2)
# plt.axis([12, 22, 12, 22])
# plt.grid(alpha=0.5)
# plt.show()

# mask = cat['gmag']-cat['rmag']<1.2
# plt.figure(figsize=(8, 8))
# plt.plot(np.arange(100), np.arange(100), color='C3')
# plt.plot(dmag_synth[mask], cat['nb_dmag'][mask], '.', ms=2)
# plt.axis([12, 22, 12, 22])
# plt.grid(alpha=0.5)
# plt.show()

cflux = cat['wiroc_flux_nb_c'] - 3/np.sqrt(cat['wiroc_flux_ivar_nb_c'])  # Define a flux less 3*sigma
cmag = 22.5 - 2.5*np.log10(cflux)  # Define a corresponding mag
lae_sel_c = (-2.25 < cmag - cmag_synth) & (cmag - cmag_synth < -0.5)  # 0.5 to 2.25 mag brighter in C than the synthetic C
print(np.sum(lae_sel_c))
lae_sel_c &= cat['wiroc_flux_nb_c'] > 3 * cat['wirod_flux_nb_d']  # At least 3x brighter in the C filter than D filter
print(np.sum(lae_sel_c))

dflux = cat['wirod_flux_nb_d'] - 3/np.sqrt(cat['wirod_flux_ivar_nb_d'])  # Define a flux less 3*sigma
dmag = 22.5 - 2.5*np.log10(dflux)  # Define a corresponding mag
lae_sel_d = (-2.25 < dmag - dmag_synth) & (dmag - dmag_synth < -0.5)  # 0.5 to 2.25 mag brighter in C than the synthetic C
print(np.sum(lae_sel_d))
lae_sel_d &= cat['wirod_flux_nb_d'] > 3 * cat['wiroc_flux_nb_c']  # At least 3x brighter in the C filter than D filter
print(np.sum(lae_sel_d))

cat[lae_sel_c].write('/global/cfs/cdirs/desi/users/rongpu/xmm_lae/wiro_c_lae_targets.fits')
cat[lae_sel_d].write('/global/cfs/cdirs/desi/users/rongpu/xmm_lae/wiro_d_lae_targets.fits')
