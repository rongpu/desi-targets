# Add sweep, photo-z and stellar mass columns
# Example:
# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python add_pz_and_stellar_mass-for_shadab.py sv3target_BGS_ANY_NBMZLS.fits

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


n_processes = 32

input_dir = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/LSScats/Target4Ang/dr9v0.57.0/sv3_v1/'
output_dir = '/global/cscratch1/sd/rongpu/target/catalogs/dr9.0/0.57.0/sv3_v1_shadab/'
stellar_mass_dir = '/global/cfs/cdirs/desi/users/rongpu/ls_dr9.0_photoz/stellar_mass'

# filelist = glob.glob('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/LSScats/Target4Ang/dr9v0.57.0/sv3/sv3target_*.fits')
# index = 0
# while index < len(filelist):
#     fn = filelist[index]
#     if "_QSO_" in fn:
#         filelist.pop(index)
#         print(fn)
#     else:
#         index+=1
# filelist = ['sv3target_BGS_ANY_NBMZLS.fits',
#              'sv3target_LRG_DES.fits',
#              'sv3target_BGS_ANY_SDECALS.fits',
#              'sv3target_BGS_ANY_NDECALS.fits',
#              'sv3target_LRG_SDECALS.fits',
#              'sv3target_ELG_SDECALS_noDES.fits',
#              'sv3target_ELG_SDECALS.fits',
#              'sv3target_ELG_DES.fits',
#              'sv3target_BGS_ANY_DES.fits',
#              'sv3target_ELG_NDECALS.fits',
#              'sv3target_LRG_NDECALS.fits',
#              'sv3target_BGS_ANY_SDECALS_noDES.fits',
#              'sv3target_ELG_NBMZLS.fits',
#              'sv3target_LRG_SDECALS_noDES.fits',
#              'sv3target_LRG_NBMZLS.fits']

fn = str(sys.argv[1])
cat_basic_path = os.path.join(input_dir, fn)
output_path = os.path.join(output_dir, fn.replace('.fits', '_more.fits'))

if os.path.isfile(output_path):
    sys.exit('File already exist: '+output_path)

cat_basic = Table(fitsio.read(cat_basic_path, columns=['RA', 'DEC', 'TARGETID', 'PHOTSYS']))

if np.all(cat_basic['PHOTSYS']=='S'):
    field = 'south'
elif np.all(cat_basic['PHOTSYS']=='N'):
    field = 'north'
else:
    raise ValueError

# #########################################################################################
# cat_basic = cat_basic[:len(cat_basic)//50]
# #########################################################################################

sweep_dir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0'.format(field)
sweep_fn_list = np.array(sorted(glob.glob(os.path.join(sweep_dir, '*.fits'))))

sweep_radec_list = [decode_sweep_name(sweep_fn) for sweep_fn in sweep_fn_list]
mask = np.array([np.any(is_in_box(cat_basic, sweep_radec)) for sweep_radec in sweep_radec_list])
print(np.sum(mask), len(mask))
sweep_fn_list = sweep_fn_list[mask]


def get_sweep_columns(sweep_fn):

    cat = Table(fitsio.read(sweep_fn, columns=['OBJID', 'BRICKID', 'RELEASE']))
    targetid = encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE'])
    idx = np.where(np.in1d(targetid, cat_basic['TARGETID']))[0]
    if len(idx)==0:
        return None
    targetid = targetid[idx]
    pz_fn = sweep_fn.replace('sweep/9.0/', 'sweep/9.0-photo-z/').replace('.fits', '-pz.fits')
    cat = Table(fitsio.read(pz_fn, rows=idx))
    cat['TARGETID'] = targetid
    cat.remove_columns(['OBJID', 'BRICKID', 'RELEASE'])

    # Add stellar mass
    stellar_mass_path = os.path.join(stellar_mass_dir, field, os.path.basename(sweep_fn).replace('.fits', '_stellar_mass.npy'))
    cat['stellar_mass'] = np.load(stellar_mass_path)[idx]

    return cat


if __name__ == '__main__':

    print('Start!')
    time_start = time.time()

    # start multiple worker processes
    with Pool(processes=n_processes) as pool:
        res = pool.map(get_sweep_columns, np.unique(sweep_fn_list))

    # Remove None elements from the list
    for index in range(len(res)-1, -1, -1):
        if res[index] is None:
            res.pop(index)

    cat_more = vstack(res, join_type='exact')
    if len(cat_more)!=len(cat_basic):
        print(len(cat_more), len(cat_basic))
        raise ValueError('different catalog length')

    # Here matching cat_more to cat_basic
    t1_reverse_sort = np.array(cat_basic['TARGETID']).argsort().argsort()
    cat_more = cat_more[np.argsort(cat_more['TARGETID'])[t1_reverse_sort]]
    if not np.all(cat_more['TARGETID']==cat_basic['TARGETID']):
        raise ValueError('different targetid')
    cat_more.remove_column('TARGETID')

    cat_more.write(output_path)

    print(time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
