from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

# cat = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/users/dstn/ODIN/xmm-N419/tractor-xmm-N419-hsc-forced.fits'))
cat = Table(fitsio.read('/Users/rongpu/Downloads/tmp/tractor-xmm-N419-hsc-forced.fits'))
print(len(cat))

mask = cat['type']!='DUP'
print('DUP', np.sum(mask)/len(mask), np.sum(mask), np.sum(~mask))
cat = cat[mask]

mask = (cat['flux_n419']>0) & (cat['flux_ivar_n419']>0)
print('Positive flux and flux_ivar in N419', np.sum(mask)/len(mask), np.sum(mask), np.sum(~mask))
cat = cat[mask]

# Extinction correction coefficients for the ODIN NB filters from Arjun's notebook
a419 = 4.3238
a501 = 3.54013
a673 = 2.43846

# Extinction-corrected fluxes and errors
cat['flux_n419_ec'] = cat['flux_n419']*10**(0.4*a419*cat['ebv'])
cat['forced_flux_g_ec'] = cat['forced_flux_g']*10**(0.4*3.240*cat['ebv'])
cat['forced_flux_r_ec'] = cat['forced_flux_r']*10**(0.4*2.276*cat['ebv'])
cat['forced_flux_i_ec'] = cat['forced_flux_i']*10**(0.4*1.633*cat['ebv'])
cat['forced_flux_z_ec'] = cat['forced_flux_z']*10**(0.4*1.263*cat['ebv'])
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    cat['flux_n419_err_ec'] = 1/np.sqrt(cat['flux_ivar_n419'])*10**(0.4*a419*cat['ebv'])
    cat['forced_flux_g_err_ec'] = 1/np.sqrt(cat['forced_flux_ivar_g'])*10**(0.4*3.240*cat['ebv'])
    cat['forced_flux_r_err_ec'] = 1/np.sqrt(cat['forced_flux_ivar_r'])*10**(0.4*2.276*cat['ebv'])
    cat['forced_flux_i_err_ec'] = 1/np.sqrt(cat['forced_flux_ivar_i'])*10**(0.4*1.633*cat['ebv'])
    cat['forced_flux_z_err_ec'] = 1/np.sqrt(cat['forced_flux_ivar_z'])*10**(0.4*1.263*cat['ebv'])

# Extinction-corrected magnitudes
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    cat['n419'] = 22.5 - 2.5*np.log10(np.clip(cat['flux_n419']*10**(0.4*a419*cat['ebv']), 1e-7, None))
    cat['gmag'] = 22.5 - 2.5*np.log10(np.clip(cat['forced_flux_g']*10**(0.4*3.240*cat['ebv']), 1e-7, None))
    cat['rmag'] = 22.5 - 2.5*np.log10(np.clip(cat['forced_flux_r']*10**(0.4*2.276*cat['ebv']), 1e-7, None))
    cat['imag'] = 22.5 - 2.5*np.log10(np.clip(cat['forced_flux_i']*10**(0.4*1.633*cat['ebv']), 1e-7, None))
    cat['zmag'] = 22.5 - 2.5*np.log10(np.clip(cat['forced_flux_z']*10**(0.4*1.263*cat['ebv']), 1e-7, None))

from astropy.coordinates import SkyCoord
from astropy import units as u

ra, dec = 35.71, -4.75

stats_search_radius = 1.62
xmm_area = np.pi * stats_search_radius**2

sky1 = SkyCoord(ra*u.degree, dec*u.degree, frame='icrs')
sky2 = SkyCoord(cat['ra']*u.degree, cat['dec']*u.degree, frame='icrs')
mask = np.array(sky2.separation(sky1).to(u.degree))<stats_search_radius
cat['in_xmm'] = mask.copy()
print(np.sum(cat['in_xmm']))

########################### LAE selection ###########################
flux_ratio = cat['forced_flux_g_ec'] / cat['flux_n419_ec']
flux_ratio_err = np.sqrt(cat['forced_flux_g_err_ec']**2 / cat['flux_n419_ec']**2 + 
                 (cat['forced_flux_g_ec'] / cat['flux_n419_ec'])**2 * (cat['flux_n419_err_ec'] / cat['flux_n419_ec'])**2)

lae_sigma_thresh = 5.
lae_n419g_thresh = 0.  # Threshold on (N419 - g) color
lae_ratio_thresh = 10**(lae_n419g_thresh/2.5)  # Threshold on FLUX_G/FLUX_N419

lae_sel = flux_ratio < lae_ratio_thresh - lae_sigma_thresh * flux_ratio_err
lae_sel &= cat['n419']<25.0
print(np.sum(lae_sel))

mask = cat['in_xmm'].copy()
mask &= lae_sel

print('LAE:')
print('Targets in DESI-XMM:', np.sum(mask))
print('Target density: {:.2f}'.format(np.sum(mask)/xmm_area))
mask &= cat['maskbits']==0
print('Target density (clean): {:.2f}'.format(np.sum(mask)/xmm_area))

mask = lae_sel.copy()
# mask &= cat['in_xmm']
mask &= cat['maskbits']==0

tt = cat[mask].copy()
print(len(tt))
tt.write('/Users/rongpu/Downloads/odin_xmm_n419_lae_targets.fits', overwrite=True)
# tt.write('/global/cfs/cdirs/desi/users/rongpu/xmm_lae/odin_xmm_n419_lae_targets.fits', overwrite=True)

