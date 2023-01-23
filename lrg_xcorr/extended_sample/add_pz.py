from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

from multiprocessing import Pool

from desitarget.targets import decode_targetid, encode_targetid


# Snippets taken from desitarget

def decode_sweep_name(sweepname):
    sweepname = os.path.basename(sweepname)

    ramin, ramax = float(sweepname[6:9]), float(sweepname[14:17])
    decmin, decmax = float(sweepname[10:13]), float(sweepname[18:21])

    if sweepname[9] == 'm':
        decmin *= -1
    if sweepname[17] == 'm':
        decmax *= -1

    return [ramin, ramax, decmin, decmax]


def is_in_box(objs, radecbox, ra_col='RA', dec_col='DEC'):

    ramin, ramax, decmin, decmax = radecbox

    # ADM check for some common mistakes.
    if decmin < -90. or decmax > 90. or decmax <= decmin or ramax <= ramin:
        msg = "Strange input: [ramin, ramax, decmin, decmax] = {}".format(radecbox)
        raise ValueError(msg)

    ii = ((objs[ra_col] >= ramin) & (objs[ra_col] < ramax)
          & (objs[dec_col] >= decmin) & (objs[dec_col] < decmax))

    return ii


def get_sweep_columns(sweep_fn, field):

    cat = Table(fitsio.read(sweep_fn, columns=['OBJID', 'BRICKID', 'RELEASE']))
    targetid = encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE'])
    if field=='north':
        idx = np.where(np.in1d(targetid, cat_basic_north['TARGETID']))[0]
    else:
        idx = np.where(np.in1d(targetid, cat_basic_south['TARGETID']))[0]
    if len(idx)==0:
        return None
    targetid = targetid[idx]

    pz_fn = os.path.basename(sweep_fn).replace('.fits', '-pz.fits')
    pz_fn = os.path.join(pz_dir, field, pz_fn)
    pz = Table(fitsio.read(pz_fn, rows=idx))
    pz['TARGETID'] = targetid
    pz.remove_columns(['OBJID', 'BRICKID', 'RELEASE'])

    return pz


print('Start!')
time_start = time.time()

pz_dir = '/global/cfs/cdirs/desi/users/rongpu/data/ls_dr9.0_desi_photoz/pz/'

pz_columns = ['Z_PHOT_MEAN', 'Z_PHOT_MEDIAN', 'Z_PHOT_STD', 'Z_PHOT_L68', 'Z_PHOT_U68', 'Z_PHOT_L95', 'Z_PHOT_U95', 'Z_SPEC', 'SURVEY', 'TRAINING']

cat_basic_path = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/dr9_extended_lrg_0.49.0_basic.fits'
cat_pz_path = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/catalogs/dr9_extended_lrg_0.49.0_pz_new.fits'

if os.path.isfile(cat_pz_path):
    sys.exit('File already exist: '+cat_pz_path)

cat_basic = Table(fitsio.read(cat_basic_path, columns=['RA', 'DEC', 'TARGETID', 'PHOTSYS']))
cat_basic_north = cat_basic[cat_basic['PHOTSYS']=='N']
cat_basic_south = cat_basic[cat_basic['PHOTSYS']=='S']

sweep_dir_north = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/north/sweep/9.0'
sweep_fn_list_north = np.array(sorted(glob.glob(os.path.join(sweep_dir_north, '*.fits'))))
sweep_radec_list_north = [decode_sweep_name(sweep_fn) for sweep_fn in sweep_fn_list_north]
mask = np.array([np.any(is_in_box(cat_basic_north, sweep_radec)) for sweep_radec in sweep_radec_list_north])
sweep_fn_list_north = np.unique(sweep_fn_list_north[mask])

sweep_dir_south = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0'
sweep_fn_list_south = np.array(sorted(glob.glob(os.path.join(sweep_dir_south, '*.fits'))))
sweep_radec_list_south = [decode_sweep_name(sweep_fn) for sweep_fn in sweep_fn_list_south]
mask = np.array([np.any(is_in_box(cat_basic_south, sweep_radec)) for sweep_radec in sweep_radec_list_south])
sweep_fn_list_south = np.unique(sweep_fn_list_south[mask])

zipped_arg_list = list(zip(sweep_fn_list_north, ['north']*len(sweep_fn_list_north)))
zipped_arg_list += list(zip(sweep_fn_list_south, ['south']*len(sweep_fn_list_south)))

n_processes = 32
with Pool(processes=n_processes) as pool:
    res = pool.starmap(get_sweep_columns, zipped_arg_list, chunksize=1)

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

cat = vstack(res, join_type='exact')
if len(cat)!=len(cat_basic):
    print(len(cat), len(cat_basic))
    raise ValueError('different catalog length')

# Here matching cat to cat_basic
t1_reverse_sort = np.array(cat_basic['TARGETID']).argsort().argsort()
cat = cat[np.argsort(cat['TARGETID'])[t1_reverse_sort]]
if not np.all(cat['TARGETID']==cat_basic['TARGETID']):
    raise ValueError('different targetid')
cat.remove_column('TARGETID')

cat = cat[pz_columns]

cat.write(cat_pz_path)

print(time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
