# Add sweep and sweep-extra columns and DR9 photo-z's
# Example:
# salloc -N 1 -C cpu -q interactive -t 4:00:00
# python add_sweep_columns_to_target_catalogs.py LRG

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
    cat = Table(fitsio.read(sweep_fn, rows=idx, columns=sweep_columns_all))
    cat['TARGETID'] = targetid

    sweep_extra_fn = sweep_fn.replace('/sweep/9.0/', '/sweep/9.0-extra/').replace('.fits', '-ex.fits')
    cat_extra = Table(fitsio.read(sweep_extra_fn, rows=idx, columns=sweep_extra_columns))

    pz_fn = sweep_fn.replace('sweep/9.0/', 'sweep/9.0-photo-z/').replace('.fits', '-pz_dr9.fits')
    pz = Table(fitsio.read(pz_fn, rows=idx))
    pz.remove_columns(['OBJID', 'BRICKID', 'RELEASE'])

    cat = hstack([cat, cat_extra, pz], join_type='exact')

    return cat


print('Start!')
time_start = time.time()

output_dir = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.1.1/resolve'

sweep_columns_1 = ['GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG', 'GAIA_ASTROMETRIC_EXCESS_NOISE', 'FITBITS',
                   'FRACFLUX_G', 'FRACFLUX_R', 'FRACFLUX_Z', 'FRACFLUX_W1', 'FRACFLUX_W2', 'FRACMASKED_G', 'FRACMASKED_R',
                   'FRACMASKED_Z', 'FRACIN_G', 'FRACIN_R', 'FRACIN_Z', 'FIBERTOTFLUX_G',
                   'SHAPE_R', 'SHAPE_R_IVAR', 'SHAPE_E1', 'SHAPE_E2', 'SERSIC', 'DCHISQ']

sweep_columns_2 = ['GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z',
                   'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'PSFDEPTH_W1', 'PSFDEPTH_W2',
                   'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z']

sweep_extra_columns = ['NEA_G', 'NEA_R', 'NEA_Z', 'BLOB_NEA_G', 'BLOB_NEA_R', 'BLOB_NEA_Z']

pz_columns = ['Z_PHOT_MEAN', 'Z_PHOT_MEDIAN', 'Z_PHOT_STD', 'Z_PHOT_L68', 'Z_PHOT_U68', 'Z_PHOT_L95', 'Z_PHOT_U95', 'Z_SPEC', 'SURVEY', 'TRAINING']

sweep_columns_all = sweep_columns_1 + sweep_columns_2
sweep_columns_all = list(set(sweep_columns_all))  # unique columns

data_dir = '/global/cfs/cdirs/desi/users/rongpu/targets/dr9.0/1.1.1/resolve'

# target_class: "LRG", "ELG", "QSO" or "BGS_ANY"
target_class = str(sys.argv[1])
target_class = target_class.upper()

print(target_class)

cat_basic_path = os.path.join(output_dir, 'dr9_{}_1.1.1_basic.fits'.format(target_class.lower()))
cat_sweep_1_path = os.path.join(output_dir, 'dr9_{}_1.1.1_sweep_1.fits'.format(target_class.lower()))
cat_sweep_2_path = os.path.join(output_dir, 'dr9_{}_1.1.1_sweep_2.fits'.format(target_class.lower()))
cat_sweep_extra_path = os.path.join(output_dir, 'dr9_{}_1.1.1_sweep_extra_1.fits'.format(target_class.lower()))
cat_pz_path = os.path.join(output_dir, 'dr9_{}_1.1.1_pz_dr9.fits'.format(target_class.lower()))

if os.path.isfile(cat_sweep_1_path):
    sys.exit('File already exist: '+cat_sweep_1_path)

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

n_processes = 128
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

cat_1 = cat[sweep_columns_1].copy()
cat_2 = cat[sweep_columns_2].copy()
cat_extra = cat[sweep_extra_columns].copy()
cat_pz = cat[pz_columns].copy()

cat_1.write(cat_sweep_1_path)
cat_2.write(cat_sweep_2_path)
cat_extra.write(cat_sweep_extra_path)
cat_pz.write(cat_pz_path)

print(time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
