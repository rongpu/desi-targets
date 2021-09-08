# Pre-compute the density of randoms around GAIA stars

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy import stats

# sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
# from match_coord import search_around, scatter_plot, match_coord


def get_density(d_ra, d_dec, d2d, plot_radius, nbins=101, min_count=None):

    bins = np.linspace(-plot_radius, plot_radius, nbins)
    bin_spacing = bins[1] - bins[0]
    bincenter = (bins[1:]+bins[:-1])/2
    mesh_ra, mesh_dec = np.meshgrid(bincenter, bincenter)
    mesh_d2d = np.sqrt(mesh_ra**2 + mesh_dec**2)
    mask = (d2d>0.01)
    count = np.histogram2d(d_ra[mask], d_dec[mask], bins=bins)[0]
    if min_count is not None:
        count[count<min_count] = np.nan
    density = count/(bin_spacing**2)
    mask = mesh_d2d >= bins.max()-bin_spacing
    density[mask] = np.nan

    return bins, density, count


randoms_dir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve'
n_randoms_catalogs = 4

gaia_path = '/global/cfs/cdirs/desi/users/rongpu/useful/gaia_edr3_g_18_dr9.fits'

# randoms_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'PHOTSYS']
gaia_columns = ['RA', 'DEC', 'PHOT_G_MEAN_MAG', 'PHOT_G_MEAN_FLUX_OVER_ERROR']

randoms_paths = sorted(glob.glob('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-[0-9]*.fits'))
randoms_paths = randoms_paths[:n_randoms_catalogs]

randoms_density = fitsio.read_header(randoms_paths[0], ext=1)['DENSITY']  # randoms per sq. deg.

##################################################################################################################################

# field = 'north'
field = str(sys.argv[1])
field = field.lower()

min_nobs = 1
maskbits = [1, 12, 13]

# unWISE maskbits: all except the SPIKE and HALO bits
wise_maskbits = [0, 2, 3, 4, 6, 7]

if field=='south':
    photsys = 'S'
else:
    photsys = 'N'

gaia = Table(fitsio.read(gaia_path, columns=gaia_columns))
print(len(gaia))

if field=='south':
    mask = (gaia['DEC']<34)
    gaia = gaia[mask]
else:
    mask = (gaia['DEC']>30)
    mask &= (gaia['RA']<310) & (gaia['RA']>75)
    gaia = gaia[mask]
print(len(gaia))

randoms = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/tmp/randoms_{}_minobs_{}_maskbits_{}.fits'.format(field, min_nobs, ''.join([str(tmp) for tmp in maskbits]))))
print('randoms:', len(randoms))

# Apply the unWISE maskbits
mask_clean = np.ones(len(randoms), dtype=bool)
for bit in wise_maskbits:
    mask_clean &= (randoms['WISEMASK_W1'] & 2**bit)==0
print('unWISE:', np.sum(~mask_clean), np.sum(~mask_clean)/len(mask_clean))
randoms = randoms[mask_clean]
print(len(randoms))

##################################################################################################################################

ra2_rand = np.array(randoms['RA'])
dec2_rand = np.array(randoms['DEC'])
sky2_rand = SkyCoord(ra2_rand*u.degree, dec2_rand*u.degree, frame='icrs')

gaia_min_list = [-np.inf, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
gaia_max_list = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]

data = {}

for index in range(len(gaia_min_list)):

    gaia_min, gaia_max = gaia_min_list[index], gaia_max_list[index]

    if index<2:
        nbins = 101
    else:
        nbins = 101

    plot_radius = 4.1
    if gaia_min==-np.inf:
        search_radius = 4000.
    else:
        search_radius = plot_radius * 1630 * 1.396**(-gaia_min)

    mask = (gaia['PHOT_G_MEAN_MAG']>gaia_min) & (gaia['PHOT_G_MEAN_MAG']<gaia_max)
    ra1 = gaia['RA'][mask]
    dec1 = gaia['DEC'][mask]
    gaia1 = gaia[mask]
    if gaia_min==-np.inf:
        title = 'GAIA_G < {:.1f}'.format(gaia_max, np.sum(mask))
    else:
        title = '{:.1f} < GAIA_G < {:.1f}'.format(gaia_min, gaia_max, np.sum(mask))

    print(title, '{} stars'.format(np.sum(mask)))

    sky1 = SkyCoord(ra1*u.degree, dec1*u.degree, frame='icrs')

    # randoms
    idx1_rand, idx2_rand, d2d_rand, _ = sky2_rand.search_around_sky(sky1, seplimit=search_radius*u.arcsec)
    print('%d nearby objects'%len(idx1_rand))
    # convert distances to numpy array in arcsec
    d2d_rand = np.array(d2d_rand.to(u.arcsec))
    gaia_g = gaia1['PHOT_G_MEAN_MAG'][idx1_rand]
    # mask_radius = 150. * 2.5**((11. - gaia_g)/3.) * 0.262
    mask_radius = 1630 * 1.396**(-gaia_g)
    d2d_rand = d2d_rand / mask_radius
    mask = d2d_rand < plot_radius
    idx1_rand = idx1_rand[mask]
    idx2_rand = idx2_rand[mask]
    d2d_rand = d2d_rand[mask]
    mask_radius = mask_radius[mask]
    d_ra_rand = (ra2_rand[idx2_rand]-ra1[idx1_rand])*3600.     # in arcsec
    d_dec_rand = (dec2_rand[idx2_rand]-dec1[idx1_rand])*3600.  # in arcsec
    # Convert d_ra_rand to actual arcsecs
    mask = d_ra_rand > 180*3600
    d_ra_rand[mask] = d_ra_rand[mask] - 360.*3600
    mask = d_ra_rand < -180*3600
    d_ra_rand[mask] = d_ra_rand[mask] + 360.*3600
    d_ra_rand = d_ra_rand * np.cos(dec1[idx1_rand]/180*np.pi)
    d_ra_rand = d_ra_rand / mask_radius
    d_dec_rand = d_dec_rand / mask_radius

    bins, density_rand, count_rand = get_density(d_ra_rand, d_dec_rand, d2d_rand, plot_radius, nbins=nbins, min_count=None)
    key_str = '{:g}_{:g}'.format(gaia_min, gaia_max)
    data[key_str+'_bins'] = bins
    data[key_str+'_density_rand'] = density_rand
    data[key_str+'_count'] = count_rand

data['n_randoms'] = len(randoms)

save_path = '/global/u2/r/rongpu/notebooks/desi_mask/data/new_mask/density_rand_gaia_{}_minobs_{}_maskbits_{}.npy'.format(field, min_nobs, ''.join([str(tmp) for tmp in maskbits]))
np.save(save_path, data)
