# find tiles -type f ! -name "*.log" > /global/u2/r/rongpu/temp/fugu_checks/everest_tiles.txt
# find exposures -type f ! -name "*.log" > /global/u2/r/rongpu/temp/fugu_checks/everest_exposures.txt
# find healpix -type f ! -name "*.log" > /global/u2/r/rongpu/temp/fugu_checks/everest_healpix.txt
# find preproc -type f ! -name "*.log" > /global/u2/r/rongpu/temp/fugu_checks/everest_preproc.txt

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

fuji_dir = '/global/cfs/cdirs/desi/spectro/redux/fuji'
guad_dir = '/global/cfs/cdirs/desi/spectro/redux/guadalupe'

############################### tiles ###############################

tiles_fn = '/global/u2/r/rongpu/temp/fugu_checks/everest_tiles.txt'
with open(tiles_fn, "r") as f:
    fn_all = f.read().splitlines()

missing_dirs = []

output_fn = '/global/u2/r/rongpu/temp/fugu_checks/tiles_missing.txt'
f = open(output_fn, 'w')

for index, fn in enumerate(fn_all):
    if index%1000==0:
        print(index, len(fn_all), '{:.1f}%'.format(100*index/len(fn_all)))
    fuji_fn = os.path.join(fuji_dir, fn)
    guad_fn = os.path.join(guad_dir, fn)
    if not (os.path.isfile(fuji_fn) or os.path.isfile(guad_fn)):
        dn = os.path.dirname(fn)
        if dn in missing_dirs:
            continue
        fuji_fn_dir = os.path.join(fuji_dir, dn)
        guad_fn_dir = os.path.join(guad_dir, dn)
        if not (os.path.isdir(fuji_fn_dir) or os.path.isdir(guad_fn_dir)):
            missing_dirs.append(dn)
            print(dn)
            f.write(dn+'\n')
        else:
            print(fn)
            f.write(fn+'\n')
f.close()

############################### exposures ###############################

exposures_fn = '/global/u2/r/rongpu/temp/fugu_checks/everest_exposures.txt'
with open(exposures_fn, "r") as f:
    fn_all = f.read().splitlines()

missing_dirs = []

output_fn = '/global/u2/r/rongpu/temp/fugu_checks/exposures_missing.txt'
f = open(output_fn, 'w')

for index, fn in enumerate(fn_all):
    if index%1000==0:
        print(index, len(fn_all), '{:.1f}%'.format(100*index/len(fn_all)))
    fuji_fn = os.path.join(fuji_dir, fn)
    guad_fn = os.path.join(guad_dir, fn)
    if not (os.path.isfile(fuji_fn) or os.path.isfile(guad_fn)):
        dn = os.path.dirname(fn)
        if dn in missing_dirs:
            continue
        fuji_fn_dir = os.path.join(fuji_dir, dn)
        guad_fn_dir = os.path.join(guad_dir, dn)
        if not (os.path.isdir(fuji_fn_dir) or os.path.isdir(guad_fn_dir)):
            missing_dirs.append(dn)
            print(dn)
            f.write(dn+'\n')
        else:
            print(fn)
            f.write(fn+'\n')
f.close()

############################### healpix ###############################

healpix_fn = '/global/u2/r/rongpu/temp/fugu_checks/everest_healpix.txt'
with open(healpix_fn, "r") as f:
    fn_all = f.read().splitlines()

missing_dirs = []

output_fn = '/global/u2/r/rongpu/temp/fugu_checks/healpix_missing.txt'
f = open(output_fn, 'w')

for index, fn in enumerate(fn_all):
    if index%1000==0:
        print(index, len(fn_all), '{:.1f}%'.format(100*index/len(fn_all)))
    fuji_fn = os.path.join(fuji_dir, fn)
    guad_fn = os.path.join(guad_dir, fn)
    if not (os.path.isfile(fuji_fn) or os.path.isfile(guad_fn)):
        dn = os.path.dirname(fn)
        if dn in missing_dirs:
            continue
        fuji_fn_dir = os.path.join(fuji_dir, dn)
        guad_fn_dir = os.path.join(guad_dir, dn)
        if not (os.path.isdir(fuji_fn_dir) or os.path.isdir(guad_fn_dir)):
            missing_dirs.append(dn)
            f.write(dn+'\n')
        else:
            print(fn)
            f.write(fn+'\n')
f.close()

############################### preproc ###############################

preproc_fn = '/global/u2/r/rongpu/temp/fugu_checks/everest_preproc.txt'
with open(preproc_fn, "r") as f:
    fn_all = f.read().splitlines()

missing_dirs = []

output_fn = '/global/u2/r/rongpu/temp/fugu_checks/preproc_missing.txt'
f = open(output_fn, 'w')

for index, fn in enumerate(fn_all):
    if index%1000==0:
        print(index, len(fn_all), '{:.1f}%'.format(100*index/len(fn_all)))
    fuji_fn = os.path.join(fuji_dir, fn)
    guad_fn = os.path.join(guad_dir, fn)
    if not (os.path.isfile(fuji_fn) or os.path.isfile(guad_fn)):
        dn = os.path.dirname(fn)
        if dn in missing_dirs:
            continue
        fuji_fn_dir = os.path.join(fuji_dir, dn)
        guad_fn_dir = os.path.join(guad_dir, dn)
        if not (os.path.isdir(fuji_fn_dir) or os.path.isdir(guad_fn_dir)):
            missing_dirs.append(dn)
            f.write(dn+'\n')
        else:
            print(fn)
            f.write(fn+'\n')
f.close()
