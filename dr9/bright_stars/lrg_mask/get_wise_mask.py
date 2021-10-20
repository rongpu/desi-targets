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

from astropy import units as u
from astropy.coordinates import SkyCoord


time_start = time.time()

output_dir = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/dev'

fn = str(sys.argv[1])
print(fn)
fn_output = os.path.basename(fn).replace('.fits', '-wisemask.npz')
output_path = os.path.join(output_dir, fn_output)

if os.path.isfile(output_path):
    sys.exit()

wise_path = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/w1_bright-2mass-lrg_mask_v1.fits'
wise = Table(fitsio.read(wise_path))
print(len(wise))

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

w1min_list = np.arange(wise['w1ab'].min()-1, wise['w1ab'].max()+1, 0.5)
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

print('Done!', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
