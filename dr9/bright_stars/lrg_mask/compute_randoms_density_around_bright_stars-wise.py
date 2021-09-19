# Pre-compute the density of randoms around WISE stars

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy import stats


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
lrgmask_dir = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/lrgmask_v1/randoms'
n_randoms_catalogs = 4

wise_path = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/w1_bright-2mass-lrg_mask_v1.fits'

randoms_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'PHOTSYS']

randoms_paths = sorted(glob.glob('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-[0-9]*.fits'))
randoms_paths = randoms_paths[:n_randoms_catalogs]

randoms_density = fitsio.read_header(randoms_paths[0], ext=1)['DENSITY']  # randoms per sq. deg.

# field = 'north'
field = str(sys.argv[1])
field = field.lower()

lrgmask_bits = [0, 1, 2, 3, 4]
min_nobs = 1

if field=='south':
    photsys = 'S'
else:
    photsys = 'N'

wise = Table(fitsio.read(wise_path))
print(len(wise))

if field=='south':
    mask = (wise['DEC']<36)
else:
    mask = (wise['DEC']>30)
    mask &= (wise['RA']<310) & (wise['RA']>75)
wise = wise[mask]
print(len(wise))

##################################################################################################################################

# randoms_stack = []
# for randoms_path in randoms_paths:
#     print(randoms_path)
#     randoms = Table(fitsio.read(randoms_path, columns=randoms_columns))
#     lrgmask_path = os.path.join(lrgmask_dir, os.path.basename(randoms_path).replace('.fits', '-lrgmask_v1.fits'))
#     lrgmask = Table(fitsio.read(lrgmask_path))
#     randoms = hstack([randoms, lrgmask], join_type='exact')
#     randoms_stack.append(randoms)
# randoms = vstack(randoms_stack)
# del randoms_stack
# randoms.write('/global/cfs/cdirs/desi/users/rongpu/tmp/randoms_lrgmask_v1.fits')

randoms = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/tmp/randoms_lrgmask_v1.fits'))

mask = (randoms['PHOTSYS']==photsys)
randoms = randoms[mask]

mask = (randoms['NOBS_G']>=min_nobs) & (randoms['NOBS_R']>=min_nobs) & (randoms['NOBS_Z']>=min_nobs)
randoms = randoms[mask]

# Apply lrgmask
randoms_clean = np.ones(len(randoms), dtype=bool)
for bit in lrgmask_bits:
    randoms_clean &= (randoms['lrg_mask'] & 2**bit)==0
randoms = randoms[randoms_clean]

# Remove pixels near the LMC
ramin, ramax, decmin, decmax = 58, 110, -90, -56
mask_remove = (randoms['RA']>ramin) & (randoms['RA']<ramax) & (randoms['DEC']>decmin) & (randoms['DEC']<decmax)
randoms = randoms[~mask_remove]
print(len(randoms))

print('randoms:', len(randoms))

##################################################################################################################################

temp = SkyCoord(ra=wise['RA']*u.degree, dec=wise['DEC']*u.degree, frame='icrs').geocentrictrueecliptic
wise['RA1'], wise['DEC1'] = np.array(temp.lon), np.array(temp.lat)

temp = SkyCoord(ra=randoms['RA']*u.degree, dec=randoms['DEC']*u.degree, frame='icrs').geocentrictrueecliptic
ra2_rand, dec2_rand = np.array(temp.lon), np.array(temp.lat)
sky2_rand = SkyCoord(ra2_rand*u.degree, dec2_rand*u.degree, frame='icrs')

w1min_list = np.arange(1., 16., 0.5)
w1max_list = w1min_list + 0.5
w1min_list = np.insert(w1min_list, 0, -np.inf)
w1max_list = np.insert(w1max_list, 0, 1.)

nomask_search_radius_list = [240, 195, 150, 135, 120, 100, 80, 70, 60, 52.5, 45, 45]

data = {}

for index in range(len(w1min_list)):

    nbins = 75

    w1min, w1max = w1min_list[index], w1max_list[index]
    mask = (wise['w1ab']>w1min) & (wise['w1ab']<=w1max)
    if np.sum(mask)==0:
        continue
    ra1, dec1 = wise['RA1'][mask], wise['DEC1'][mask]
    sky1 = SkyCoord(ra1*u.degree, dec1*u.degree, frame='icrs')
    if index>=len(w1min_list)-len(nomask_search_radius_list):
        search_radius = nomask_search_radius_list[index-(len(w1min_list)-len(nomask_search_radius_list))]
    else:
        search_radius = wise['radius'][mask].max() * 4.1

    if w1min==-np.inf:
        title = 'WISE_W1_AB <= {:.1f}'.format(w1max, np.sum(mask))
    else:
        title = '{:.1f} < WISE_W1_AB <= {:.1f}'.format(w1min, w1max, np.sum(mask))

    print(title, '{} stars'.format(np.sum(mask)))

    # randoms
    idx1_rand, idx2_rand, d2d_rand, _ = sky2_rand.search_around_sky(sky1, seplimit=search_radius*u.arcsec)
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

save_path = '/global/u2/r/rongpu/notebooks/desi_mask/data/lrgmask_v1/density_rand_wise_{}_minobs_{}_lrgmask_{}.npy'.format(field, min_nobs, ''.join([str(tmp) for tmp in lrgmask_bits]))
np.save(save_path, data)
