from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from desitarget.targets import decode_targetid, encode_targetid

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


odin_hip = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/xmm_lae/odin_xmm_n419_lae_targets_high_priority_20221215.fits'))
odin_lop = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/xmm_lae/odin_xmm_n419_lae_targets_low_priority_20221215.fits'))
wiroc = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/xmm_lae/wiro_c_lae_targets.fits'))
wirod = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/xmm_lae/wiro_d_lae_targets.fits'))

odin_hip['targetid'] = encode_targetid(odin_hip['objid'], odin_hip['brickid'], odin_hip['release'])
odin_lop['targetid'] = encode_targetid(odin_lop['objid'], odin_lop['brickid'], odin_lop['release'])
wiroc['targetid'] = encode_targetid(wiroc['objid'], wiroc['brickid'], wiroc['release'])
wirod['targetid'] = encode_targetid(wirod['objid'], wirod['brickid'], wirod['release'])

print(len(odin_hip), len(odin_hip)-len(np.unique(odin_hip['targetid'])))
print(len(odin_lop), len(odin_lop)-len(np.unique(odin_lop['targetid'])))
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

odin = vstack([odin_hip, odin_lop])
odin['odin_hip'] = np.in1d(odin['targetid'], odin_hip['targetid'])
odin['odin_lop'] = np.in1d(odin['targetid'], odin_lop['targetid'])
print(len(odin), len(np.unique(odin['targetid'])))
if len(odin)!=len(np.unique(odin['targetid'])):
    _, idx = np.unique(odin['targetid'], return_index=True)
    idx = np.sort(idx)
    odin = odin[idx]
    print(len(odin), len(np.unique(odin['targetid'])))

# Assign goal exposure time (in number of 15-minute passes)
odin['n_pass'] = 0
odin['n_pass'][odin['n419']<24.25] = 4
odin['n_pass'][odin['n419']>=24.25] = 8
wiro['n_pass'] = 8

# Remove duplicates
wiro['odin_hip'] = False
wiro['odin_lop'] = False
idx1, idx2, d2d, d_ra, d_dec = match_coord(wiro['ra'], wiro['dec'], odin['ra'], odin['dec'], search_radius=1., plot_q=True)
wiro['odin_hip'][idx1] = odin['odin_hip'][idx2]
wiro['odin_lop'][idx1] = odin['odin_lop'][idx2]
mask = np.in1d(np.arange(len(odin)), idx2)
odin = odin[~mask]
print(len(odin))

wiro_columns = ['ra', 'dec', 'wiro_c', 'wiro_d', 'odin_hip', 'odin_lop', 'n_pass']
odin_columns = ['ra', 'dec', 'odin_hip', 'odin_lop', 'n_pass']
cat = vstack([wiro[wiro_columns], odin[odin_columns]]).filled(False)

cat['priority'] = '      '
mask = cat['wiro_c'] | cat['wiro_d']
cat['priority'][mask] = 'high'
mask = cat['odin_hip'].copy()
cat['priority'][mask] = 'medium'
mask = cat['odin_lop'].copy()
cat['priority'][mask] = 'low'

cat.write('/global/cfs/cdirs/desi/users/rongpu/xmm_lae/xmm_odin_wiro_merged_targets.fits')

