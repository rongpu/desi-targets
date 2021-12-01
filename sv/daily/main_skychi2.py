# Get the {G,R,Z}SKYCHI2PDF values from exposure-qa/PETALQA
# for the investigation of LRG redshift efficiency vs speed

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits


coadd_type = 'cumulative'

tiles = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
print(len(tiles), len(np.unique(tiles['TILEID'])))

mask = tiles['SURVEY']=='main'
mask &= tiles['FAPRGRM']=='dark'
tiles = tiles[mask]
print(len(tiles))

# mask = tiles['OBSSTATUS']=='obsend'
mask = tiles['QA']=='good'
tiles = tiles[mask]
print(len(tiles))

mask = tiles['LASTNIGHT']>20210801  # post shutdown
tiles = tiles[mask]
print(len(tiles))

columns_3 = ['TILEID', 'NIGHT', 'EXPID']

tileid_list = tiles['TILEID']

data_dir = '/global/cfs/cdirs/desi/spectro/redux/daily/tiles/{}'.format(coadd_type)

cat_stack = []

for tileid in tileid_list:

    print(tileid)

    fn_list = sorted(glob.glob(os.path.join(data_dir, str(tileid), '*/redrock-*.fits')))
    # for fn in fn_list:
    fn = fn_list[0]
    tmp3 = Table(fitsio.read(fn, ext=3, columns=columns_3))
    _, idx = np.unique(tmp3['EXPID'], return_index=True)
    tmp3 = tmp3[idx]
    cat_stack.append(tmp3)

cat = vstack(cat_stack)
print()
print(len(cat))

qa_all = []

for index in range(len(cat)):
    tileid = cat['TILEID'][index]
    night = cat['NIGHT'][index]
    expid = cat['EXPID'][index]
    expid_str = str(expid).zfill(8)
    fn = '/global/cfs/cdirs/desi/spectro/redux/daily/exposures/{0}/{1}/exposure-qa-{1}.fits'.format(night, expid_str)
    # print(os.path.isfile(fn))
    qa = Table(fitsio.read(fn, ext='PETALQA'))
    qa['TILEID'] = tileid
    qa['NIGHT'] = night
    qa['EXPID'] = expid
    qa_all.append(qa)
qa = vstack(qa_all)

qa.write('/global/cfs/cdirs/desi/users/rongpu/spectro/daily/petal_qa_{}_20211028.fits'.format(coadd_type))

