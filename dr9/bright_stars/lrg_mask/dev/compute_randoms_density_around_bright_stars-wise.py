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

from scipy.interpolate import interp1d

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

wise_path = '/global/cfs/cdirs/desi/users/rongpu/useful/w1_bright-2mass-13.3-dr9.fits'

randoms_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'PHOTSYS', 'WISEMASK_W1']

randoms_paths = sorted(glob.glob('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-[0-9]*.fits'))
randoms_paths = randoms_paths[:n_randoms_catalogs]

randoms_density = fitsio.read_header(randoms_paths[0], ext=1)['DENSITY']  # randoms per sq. deg.

##################################################################################################################################

# # v1
# w1_mags = [0, 0.5, 1, 2, 2.5, 4.0, 5.0, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0]
# w1_radii = [600, 600, 500, 425, 375, 350, 350, 300, 250, 200, 175, 150, 100, 100]
# w1_max_mag = 9.0

# v2
w1_mags = [0, 0.5, 1, 1.5, 2, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
w1_radii = [600, 600, 550, 500, 475, 425, 400, 400, 375, 375, 375, 350, 325, 300, 250, 225, 175, 150, 75, 50, 50]
w1_max_mag = 10.0

f_radius = interp1d(w1_mags, w1_radii, bounds_error=False, fill_value='extrapolate')

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

wise = Table(fitsio.read(wise_path))
print(len(wise))

if field=='south':
    mask = (wise['DEC']<36)
else:
    mask = (wise['DEC']>31)
    mask &= (wise['RA']<310) & (wise['RA']>75)
wise = wise[mask]
print(len(wise))

wise['w1ab'] = np.array(wise['W1MPRO']) + 2.699

wise['radius'] = f_radius(wise['w1ab'])

# randoms_stack = []

# for randoms_path in randoms_paths:
#     print(randoms_path)

#     randoms = Table(fitsio.read(randoms_path, columns=randoms_columns))

#     mask = (randoms['PHOTSYS']==photsys)
#     randoms = randoms[mask]

#     mask = (randoms['NOBS_G']>=min_nobs) & (randoms['NOBS_R']>=min_nobs) & (randoms['NOBS_Z']>=min_nobs)
#     randoms = randoms[mask]

#     # Apply masks
#     randoms_clean = np.ones(len(randoms), dtype=bool)
#     for bit in maskbits:
#         randoms_clean &= (randoms['MASKBITS'] & 2**bit)==0
#     randoms = randoms[randoms_clean]

#     # Remove pixels near the LMC
#     ramin, ramax, decmin, decmax = 58, 110, -90, -56
#     mask_remove = (randoms['RA']>ramin) & (randoms['RA']<ramax) & (randoms['DEC']>decmin) & (randoms['DEC']<decmax)
#     randoms = randoms[~mask_remove]
#     print(len(randoms))

#     randoms_stack.append(randoms)

# randoms = vstack(randoms_stack)
# del randoms_stack
# randoms.write('/global/cfs/cdirs/desi/users/rongpu/tmp/randoms_{}_minobs_{}_maskbits_{}.fits'.format(field, min_nobs, ''.join([str(tmp) for tmp in maskbits])))

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

temp = SkyCoord(ra=wise['RA']*u.degree, dec=wise['DEC']*u.degree, frame='icrs').geocentrictrueecliptic
wise['RA1'], wise['DEC1'] = np.array(temp.lon), np.array(temp.lat)

temp = SkyCoord(ra=randoms['RA']*u.degree, dec=randoms['DEC']*u.degree, frame='icrs').geocentrictrueecliptic
ra2_rand, dec2_rand = np.array(temp.lon), np.array(temp.lat)
sky2_rand = SkyCoord(ra2_rand*u.degree, dec2_rand*u.degree, frame='icrs')

################################################## Remove the masked objects ########################################################

w1min_list = np.arange(-3, w1_max_mag+1, 0.5)
w1max_list = w1min_list + 0.5

idx2_remove = []

for index in range(len(w1min_list)):

    w1min, w1max = w1min_list[index], w1max_list[index]
    mask = (wise['w1ab']>w1min) & (wise['w1ab']<w1max)
    mask &= wise['w1ab']<w1_max_mag

    if np.sum(mask)==0:
        continue
    ra1, dec1 = wise['RA1'][mask], wise['DEC1'][mask]
    sky1 = SkyCoord(ra1*u.degree, dec1*u.degree, frame='icrs')
    wise1 = wise[mask].copy()
    search_radius = wise1['radius'].max()
    
    if w1min==-np.inf:
        title = 'WISE_W1_AB < {:.1f}'.format(w1max, np.sum(mask))
    else:
        title = '{:.1f} < WISE_W1_AB < {:.1f}'.format(w1min, w1max, np.sum(mask))

    print(title, '{} stars'.format(np.sum(mask)))

    #############################################################################

    # Objects
    idx1_rand, idx2_rand, d2d_rand, _ = sky2_rand.search_around_sky(sky1, seplimit=search_radius*u.arcsec)
    if len(idx1_rand)==0:
        continue
    print('{} nearby objects around {} stars'.format(len(np.unique(idx2_rand)), len(np.unique(idx1_rand))))
    d2d_rand = np.array(d2d_rand.to(u.arcsec))  # convert distances to numpy array in arcsec

    mask = d2d_rand < wise1['radius'][idx1_rand]
    idx2_remove.append(np.unique(idx2_rand[mask]))

idx2_remove = np.unique(np.concatenate(idx2_remove))
mask_remove = np.full(len(randoms), False)
mask_remove[idx2_remove] = True
print('New mask:', np.sum(mask_remove), np.sum(mask_remove)/len(mask_remove))

randoms['w1_mask'] = mask_remove.copy()
randoms.write('/global/cfs/cdirs/desi/users/rongpu/tmp/randoms_{}_minobs_{}_maskbits_{}-wise_v2.fits'.format(field, min_nobs, ''.join([str(tmp) for tmp in maskbits])))

##################################################################################################################################

w1min_list = np.arange(0., 12., 0.5)
w1max_list = w1min_list + 0.5
search_radius_list = 221.3 * 1.393**(-(w1max_list-2.599-2.699)) * 5.
mask = w1min_list<4
search_radius_list[mask] = search_radius_list[mask] * 3/5

data = {}

for index in range(len(w1min_list)):

    if index<2:
        nbins = 75
    else:
        nbins = 75
    search_radius = search_radius_list[index]

    w1min, w1max = w1min_list[index], w1max_list[index]
    mask = (wise['w1ab']>w1min) & (wise['w1ab']<w1max)
    ra1, dec1 = wise['RA1'][mask], wise['DEC1'][mask]
    sky1 = SkyCoord(ra1*u.degree, dec1*u.degree, frame='icrs')

    if w1min==-np.inf:
        title = 'WISE_W1_AB < {:.1f}'.format(w1max, np.sum(mask))
    else:
        title = '{:.1f} < WISE_W1_AB < {:.1f}'.format(w1min, w1max, np.sum(mask))

    print(title, '{} stars'.format(np.sum(mask)))

    # randoms
    idx1_rand, idx2_rand, d2d_rand, _ = sky2_rand.search_around_sky(sky1, seplimit=search_radius*u.arcsec)
    mask = ~np.in1d(idx2_rand, idx2_remove)
    idx1_rand = idx1_rand[mask]
    idx2_rand = idx2_rand[mask]
    d2d_rand = d2d_rand[mask]
    if len(idx1_rand)==0:
        continue
    print('{} nearby objects around {} stars'.format(len(np.unique(idx2_rand)), len(np.unique(idx1_rand))))
    d2d_rand = np.array(d2d_rand.to(u.arcsec))  # convert distances to numpy array in arcsec
    d_ra_rand = (ra2_rand[idx2_rand]-ra1[idx1_rand])*3600.     # in arcsec
    d_dec_rand = (dec2_rand[idx2_rand]-dec1[idx1_rand])*3600.  # in arcsec
    # Convert d_ra_rand to actual arcsecs
    mask = d_ra_rand > 180*3600
    d_ra_rand[mask] = d_ra_rand[mask] - 360.*3600
    mask = d_ra_rand < -180*3600
    d_ra_rand[mask] = d_ra_rand[mask] + 360.*3600
    d_ra_rand = d_ra_rand * np.cos(dec1[idx1_rand]/180*np.pi)

    bins, density_rand, count_rand = get_density(d_ra_rand, d_dec_rand, d2d_rand, search_radius, nbins=nbins, min_count=None)
    key_str = '{:g}_{:g}'.format(w1min, w1max)
    data[key_str+'_bins'] = bins
    data[key_str+'_density_rand'] = density_rand
    data[key_str+'_count'] = count_rand

data['n_randoms'] = len(randoms)

save_path = '/global/u2/r/rongpu/notebooks/desi_mask/data/new_mask/density_rand_wise_{}_minobs_{}_maskbits_{}_wise_v2.npy'.format(field, min_nobs, ''.join([str(tmp) for tmp in maskbits]))
np.save(save_path, data)
