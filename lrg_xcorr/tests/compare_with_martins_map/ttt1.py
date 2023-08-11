from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp

sys.path.append(os.path.expanduser('~/git/desi-examples/imaging_systematics'))
from plot_healpix_map import plot_map


params = {'legend.fontsize': 'x-large',
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)

plt.rcParams['image.cmap'] = 'jet'

weighted = True

min_nobs = 2
maskbits = []
custom_mask_name = 'lrgmask_v1.1'

mask_str = ''.join([str(tmp) for tmp in maskbits])
if custom_mask_name!='':
    mask_str += '_' + custom_mask_name

randoms_counts_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/resolve/counts'
randoms_systematics_dir = '/global/cfs/cdirs/desi/users/rongpu/data/imaging_sys/randoms_stats/0.49.0/resolve/systematics'
target_densities_dir = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/density_maps/main_lrg'

if weighted:
    target_densities_dir = os.path.join(target_densities_dir, 'linear_weights')

dpi_dict = {64: 100, 128: 100, 256: 300, 512: 1600}
xsize_dict = {64: 4000, 128: 4000, 256: 6000, 512: 16000}
vrange_dict = {64: 8, 128: 12, 256: 18}
vcenter_dict = {1: 81.0, 2: 147.5, 3: 164.2, 4: 149.4}  # average density in DECam imaging
# vrange_dict = {64: [0, 1200], 128: [-200, 1400], 256: [-600, 1800]}

nside = 64

for bin_index in range(1, 5):
    hp_table = Table(fitsio.read('/pscratch/sd/r/rongpu/tmp/lrg_bin_{}_{}.fits'.format(bin_index, nside)))

    vcenter = vcenter_dict[bin_index]
    vrange = int(np.sqrt(vcenter)*vrange_dict[nside])
    vmin, vmax = vcenter - vrange, vcenter + vrange
    if weighted:
        vmin, vmax = vmin/vcenter, vmax/vcenter

    plot_dir = '/pscratch/sd/r/rongpu/tmp/lrg_density_maps'

    plot_path = os.path.join(plot_dir, 'density_lrg_pz_bin_{}_{}_martin.png'.format(bin_index, nside))
    plot_map(nside, hp_table['delta']+1, pix=hp_table['HPXPIXEL'],
             vmin=vmin, vmax=vmax, xsize=xsize_dict[nside], dpi=dpi_dict[nside],
             title='LRG photo-z bin {} NSIDE={}'.format(bin_index, nside), save_path=plot_path, show=False)
