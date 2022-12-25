from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/moregit/desitarget/py/'))
sys.path.append(os.path.expanduser('~/moregit/desiutil/py/'))
from desitarget.targets import decode_targetid, encode_targetid

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


odin = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/xmm_lae/odin_xmm_n419_lae_targets.fits'))
wiroc = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/xmm_lae/wiro_c_lae_targets.fits'))
wirod = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/xmm_lae/wiro_d_lae_targets.fits'))
# odin = Table(fitsio.read('/Users/rongpu/Downloads/odin_wiro_data/odin_xmm_n419_lae_targets.fits'))
# wiroc = Table(fitsio.read('/Users/rongpu/Downloads/odin_wiro_data/wiro_c_lae_targets.fits'))
# wirod = Table(fitsio.read('/Users/rongpu/Downloads/odin_wiro_data/wiro_d_lae_targets.fits'))

odin['targetid'] = encode_targetid(odin['objid'], odin['brickid'], odin['release'])
wiroc['targetid'] = encode_targetid(wiroc['objid'], wiroc['brickid'], wiroc['release'])
wirod['targetid'] = encode_targetid(wirod['objid'], wirod['brickid'], wirod['release'])

print(len(odin), len(odin)-len(np.unique(odin['targetid'])))
print(len(wiroc), len(wiroc)-len(np.unique(wiroc['targetid'])))
print(len(wirod), len(wirod)-len(np.unique(wirod['targetid'])))

wiro = vstack([wiroc, wirod])
wiro['wiro_c'] = np.in1d(wiro['targetid'], wiroc['targetid'])
wiro['wiro_d'] = np.in1d(wiro['targetid'], wirod['targetid'])
print(len(wiro), len(np.unique(wiro['targetid'])))
if len(wiro)!=len(np.unique(wiro['targetid'])):
    _, idx = np.unique(wiro['targetid'], return_index=True)
    idx = np.sort(idx)
    wiro = wiro[idx]
    print(len(wiro), len(np.unique(wiro['targetid'])))

odin['odin'] = True
odin['odin_bright'] = odin['n419']<24.25
odin['odin_faint'] = odin['n419']>=24.25
print(len(odin), len(np.unique(odin['targetid'])))
if len(odin)!=len(np.unique(odin['targetid'])):
    _, idx = np.unique(odin['targetid'], return_index=True)
    idx = np.sort(idx)
    odin = odin[idx]
    print(len(odin), len(np.unique(odin['targetid'])))

odin['fa_class'] = '           '
odin['fa_class'][odin['n419']<24.25] = 'ODIN_BRIGHT'
odin['fa_class'][odin['n419']>=24.25] = 'ODIN_FAINT'
wiro['fa_class'] = 'WIRO'

# Remove duplicates
wiro['odin'] = False
idx1, idx2, d2d, d_ra, d_dec = match_coord(wiro['ra'], wiro['dec'], odin['ra'], odin['dec'], search_radius=1., plot_q=True)
wiro['odin'][idx1] = True
wiro['odin_bright'] = False
wiro['odin_bright'][idx1] = odin['odin_bright'][idx2]
wiro['odin_faint'] = False
wiro['odin_faint'][idx1] = odin['odin_faint'][idx2]
mask = np.in1d(np.arange(len(odin)), idx2)
odin = odin[~mask]
print(len(odin))

wiro_columns = ['ra', 'dec', 'wiro_c', 'wiro_d', 'odin', 'odin_bright', 'odin_faint', 'fa_class']
odin_columns = ['ra', 'dec', 'odin', 'odin_bright', 'odin_faint', 'fa_class']
cat = vstack([wiro[wiro_columns], odin[odin_columns]]).filled(False)

cat.write('/global/cfs/cdirs/desi/users/rongpu/xmm_lae/xmm_odin_wiro_merged_targets_1.fits', overwrite=True)
# cat.write('/Users/rongpu/Downloads/odin_wiro_data/xmm_odin_wiro_merged_targets.fits', overwrite=True)

