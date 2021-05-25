# Create cross-matched catalogs of LRGs in HSC

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

from desitarget.targets import decode_targetid, encode_targetid


for field in ['south', 'north']:

    data_dir = '/global/cfs/cdirs/desi/target/analysis/truth/dr9.0/{}/matched/'.format(field)

    filelist = ['hsc-pdr2-wide-w01-reduced-match.fits',
                'hsc-pdr2-wide-w02-reduced-match.fits',
                'hsc-pdr2-wide-w03-reduced-match.fits',
                'hsc-pdr2-wide-w04-reduced-match.fits',
                'hsc-pdr2-wide-w05-reduced-match.fits',
                'hsc-pdr2-wide-w06-reduced-match.fits',
                'hsc-pdr2-wide-w07-reduced-match.fits',
    ]
    lrg = Table(fitsio.read('/global/cscratch1/sd/rongpu/target/catalogs/dr9.0/0.49.0/dr9_sv1_lrg_{}_0.49.0_basic.fits'.format(field), columns=['TARGETID']))

    hsc_stack = []
    ls_stack = []
    for fn in filelist:
        hsc_path = os.path.join(data_dir, fn)
        ls_path = hsc_path.replace('/matched/', '/matched/ls-dr9.0-')
        if not os.path.isfile(hsc_path):
            continue
        cat = Table(fitsio.read(ls_path, columns=['OBJID', 'BRICKID', 'RELEASE']))
        targetid = encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE'])
        mask = np.in1d(targetid, lrg['TARGETID'])
        targetid = targetid[mask]
        print(np.sum(mask)/len(mask), len(cat))
        idx = np.where(mask)[0]
        if len(idx)>0:
            hsc = Table(fitsio.read(hsc_path, rows=idx))
            ls = Table(fitsio.read(ls_path, rows=idx))
            ls['TARGETID'] = targetid
            hsc_stack.append(hsc)
            ls_stack.append(ls)
    hsc = vstack(hsc_stack)
    ls = vstack(ls_stack)

    hsc.write('/global/cscratch1/sd/rongpu/target/catalogs/dr9.0/0.49.0/hsc/sv1_lrg_{}_hsc.fits'.format(field), overwrite=True)
    ls.write('/global/cscratch1/sd/rongpu/target/catalogs/dr9.0/0.49.0/hsc/sv1_lrg_{}_ls.fits'.format(field), overwrite=True)

##################################################################################################################################

