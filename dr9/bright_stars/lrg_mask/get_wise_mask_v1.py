# salloc -N 1 -C haswell -q interactive -t 4:00:00
# parallel --jobs 6 < get_wise_mask.txt; exit
# timing: 3min * 200 / 6 / 60

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

from scipy.interpolate import interp1d
from astropy import units as u
from astropy.coordinates import SkyCoord


output_dir = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/dev'

fn = str(sys.argv[1])
print(fn)
fn_output = os.path.basename(fn).replace('.fits', '-wisemask.npz')
output_path = os.path.join(output_dir, fn_output)

if os.path.isfile(output_path):
    sys.exit()

# WISE mask v4
w1_mags = [0, 0.5, 1, 1.5, 2, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
w1_radii = [600, 600, 550, 500, 475, 425, 400, 400, 390, 392.5, 395, 370, 360, 330, 275, 240, 210, 165, 100, 75, 60]
w1_max_mag = 10.0

f_radius = interp1d(w1_mags, w1_radii, bounds_error=False, fill_value='extrapolate')

wise_path = '/global/cfs/cdirs/desi/users/rongpu/useful/w1_bright-2mass-13.3-dr9.fits'
wise = Table(fitsio.read(wise_path))
# print(len(wise))

wise['w1ab'] = np.array(wise['W1MPRO']) + 2.699
mask = wise['w1ab']<w1_max_mag
wise = wise[mask]
# print(len(wise))

wise['radius'] = f_radius(wise['w1ab'])

# columns = ['RA', 'DEC', 'MASKBITS', 'WISEMASK_W1', 'PHOTSYS']
columns = ['RA', 'DEC']
hdu = fits.open(fn)
cat = Table()
for col in columns:
    cat[col] = np.copy(hdu[1].data[col])

# ############################
# cat = cat[:len(cat)//10]
# ############################

ra2, dec2 = np.array(cat['RA']), np.array(cat['DEC'])
sky2 = SkyCoord(ra2*u.degree, dec2*u.degree, frame='icrs')

del cat

#############################################################################

w1min_list = np.arange(wise['w1ab'].min()-1, w1_max_mag+1, 0.5)
w1max_list = w1min_list + 0.5

idx2_remove = []

for index in range(len(w1min_list)):

    w1min, w1max = w1min_list[index], w1max_list[index]
    mask = (wise['w1ab']>w1min) & (wise['w1ab']<w1max)

    if np.sum(mask)==0:
        continue
    ra1, dec1 = wise['RA'][mask], wise['DEC'][mask]
    sky1 = SkyCoord(ra1*u.degree, dec1*u.degree, frame='icrs')
    wise1 = wise[mask].copy()
    search_radius = wise1['radius'].max()
    
    if w1min==-np.inf:
        title = 'WISE_W1_AB < {:.1f}'.format(w1max, np.sum(mask))
    else:
        title = '{:.1f} < WISE_W1_AB < {:.1f}'.format(w1min, w1max, np.sum(mask))

    # print(title, '{} stars'.format(np.sum(mask)))

    # Objects
    idx1, idx2, d2d, _ = sky2.search_around_sky(sky1, seplimit=search_radius*u.arcsec)
    if len(idx1)==0:
        continue
    # print('{} nearby objects around {} stars'.format(len(np.unique(idx2)), len(np.unique(idx1))))
    d2d = np.array(d2d.to(u.arcsec))  # convert distances to numpy array in arcsec

    mask = d2d < wise1['radius'][idx1]
    idx2_remove.append(np.unique(idx2[mask]))

idx2_remove = np.unique(np.concatenate(idx2_remove))
mask_remove = np.full(len(ra2), False)
mask_remove[idx2_remove] = True
print('New mask:', np.sum(mask_remove), np.sum(mask_remove)/len(mask_remove))

data = {}
data['wise_mask'] = mask_remove.copy()
np.savez_compressed(output_path, **data)

