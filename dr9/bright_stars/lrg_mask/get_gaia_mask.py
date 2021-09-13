# salloc -N 1 -C haswell -q interactive -t 4:00:00
# parallel --jobs 5 < get_gaia_mask.txt; exit
# timing: 7min * 200 / 5 / 60

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
fn_output = os.path.basename(fn).replace('.fits', '-gaiamask.npz')
output_path = os.path.join(output_dir, fn_output)

if os.path.isfile(output_path):
    sys.exit()

gaia_path = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_lrg_mask_v1.fits'
gaia_columns = ['RA', 'DEC', 'mask_mag', 'radius_south', 'radius_north']
max_mag = 18.

hdu = fits.open(gaia_path)
gaia = Table()
for col in gaia_columns:
    gaia[col] = np.copy(hdu[1].data[col])

####################### Add custom masks ########################

custom_mask_fn = '/global/cfs/cdirs/desi/users/rongpu/misc/desi_custom_mask.txt'
with open(custom_mask_fn, 'r') as f:
    lines = list(map(str.strip, f.readlines()))
    
# circular mask
ra, dec, radius = [], [], []
# rectangular mask
ramin, ramax, decmin, decmax = [], [], [], []

circ_mask_data = []
rect_mask_data = []

for line in lines:
    if line!='' and line[0]!='#':
        line = line[:line.find('#')]
        line = list(map(float, line.split(',')))
        if len(line)==3:
            circ_mask_data.append(line)
        elif len(line)==4:
            rect_mask_data.append(line)
        else:
            raise ValueError

circ_mask_data = np.array(circ_mask_data)
rect_mask_data = np.array(rect_mask_data)
print('Custom circular mask:', len(circ_mask_data))
print('Custom rectangular mask:', len(rect_mask_data))

cm = Table()
cm['RA'] = circ_mask_data[:, 0]
cm['DEC'] = circ_mask_data[:, 1]
cm['mask_mag'] = np.nan
cm['radius_south'] = circ_mask_data[:, 2]
cm['radius_north'] = circ_mask_data[:, 2]

gaia = vstack([gaia, cm], join_type='exact')
print(len(gaia))

##############################################################

columns = ['RA', 'DEC', 'PHOTSYS']
hdu = fits.open(fn)
cat = Table()
for col in columns:
    cat[col] = np.copy(hdu[1].data[col])

# ############################
# cat = cat[:len(cat)//10]
# ############################

ra2, dec2 = np.array(cat['RA']), np.array(cat['DEC'])
sky2 = SkyCoord(ra2*u.degree, dec2*u.degree, frame='icrs')

#############################################################################

gaia_min_list = np.arange(np.floor(gaia['mask_mag'][np.isfinite(gaia['mask_mag'])].min())-1, max_mag+1, 0.5)
gaia_max_list = gaia_min_list + 0.5

mask_remove_north = np.full(len(ra2), False)
mask_remove_south = np.full(len(ra2), False)

mask_remove_north_bright = np.full(len(ra2), False)
mask_remove_south_bright = np.full(len(ra2), False)

for index in range(len(gaia_min_list)):

    gaia_min, gaia_max = gaia_min_list[index], gaia_max_list[index]
    mask = (gaia['mask_mag']>gaia_min) & (gaia['mask_mag']<gaia_max)
    if index==0:
        mask |= np.isnan(gaia['mask_mag'])

    if np.sum(mask)==0:
        continue
    ra1, dec1 = gaia['RA'][mask], gaia['DEC'][mask]
    sky1 = SkyCoord(ra1*u.degree, dec1*u.degree, frame='icrs')
    gaia1 = gaia[mask].copy()
    search_radius = np.maximum(gaia1['radius_north'].max(), gaia1['radius_south'].max())
    
    if gaia_min==-np.inf:
        title = 'GAIA_G < {:.1f}'.format(gaia_max, np.sum(mask))
    else:
        title = '{:.1f} < GAIA_G < {:.1f}'.format(gaia_min, gaia_max, np.sum(mask))

    # print(title, '{} stars'.format(np.sum(mask)))

    # Objects
    idx1, idx2, d2d, _ = sky2.search_around_sky(sky1, seplimit=search_radius*u.arcsec)
    if len(idx1)==0:
        continue
    # print('{} nearby objects around {} stars'.format(len(np.unique(idx2)), len(np.unique(idx1))))
    d2d = np.array(d2d.to(u.arcsec))  # convert distances to numpy array in arcsec

    mask_north = d2d < gaia1['radius_north'][idx1]
    mask_remove_north[idx2[mask_north]] = True

    mask_south = d2d < gaia1['radius_south'][idx1]
    mask_remove_south[idx2[mask_south]] = True

    mask_bright = gaia1['mask_mag'][idx1]<16
    mask_remove_north_bright[idx2[mask_north & mask_bright]] = True
    mask_remove_south_bright[idx2[mask_south & mask_bright]] = True

mask_north = cat['PHOTSYS']=='N'
mask_south = cat['PHOTSYS']=='S'

mask_remove_north &= mask_north
mask_remove_south &= mask_south
mask_remove_north_bright &= mask_north
mask_remove_south_bright &= mask_south

mask_remove = mask_remove_south | mask_remove_north
mask_remove_bright = mask_remove_south_bright | mask_remove_north_bright

####################### Rectangular masks ########################

mask_rect = np.full(len(cat), False)
for radec in rect_mask_data:
    ramin, ramax, decmin, decmax = radec
    mask_rect |= (cat['RA']>ramin) & (cat['RA']<ramax) & (cat['DEC']>decmin) & (cat['DEC']<decmax)
mask_remove |= mask_rect
mask_remove_bright |= mask_rect

###################################################################

print('New mask:', np.sum(mask_remove), np.sum(mask_remove)/len(mask_remove))
print('New mask (bright):', np.sum(mask_remove_bright), np.sum(mask_remove_bright)/len(mask_remove_bright))

data = {}
data['gaia_mask'] = mask_remove.copy()
data['gaia_bright_mask'] = mask_remove_bright.copy()
np.savez_compressed(output_path, **data)

print('Done!', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
