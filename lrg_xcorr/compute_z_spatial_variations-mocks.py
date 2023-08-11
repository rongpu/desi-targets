from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp
from multiprocessing import Pool

nmad = lambda x: 1.4826 * np.median(np.abs(x-np.median(x)))


nside = 16

for mock_subset in range(25):

    cat = Table(fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/mocks/FirstGenMocks/AbacusSummit/Y1v1/mock{}/LSScats/LRG_full.dat.fits'.format(mock_subset)))
    mask = np.isfinite(cat['Z_not4clus'])
    cat = cat[mask]
    print(len(cat))
    cat.rename_column('Z_not4clus', 'Z')

    tt = Table.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/dndz/iron_v0.2/main_lrg_pz_dndz_iron_v0.2_dz_0.01.txt', format='ascii.commented_header')

    y1_area = 3756

    np.random.seed(781273+mock_subset)
    mask_all = np.full(len(cat), False)
    mask_dict = {}
    for bin_index in range(1, 5):
        mask = np.full(len(cat), False)
        for z_index in range(len(tt)):
            zmin, zmax = tt['zmin'][z_index], tt['zmax'][z_index]
            idx = np.where(~mask_all & (cat['Z']>=zmin) & (cat['Z']<zmax))[0]
            n_samp = int(np.round(tt['bin_{}_combined'.format(bin_index)][z_index]*y1_area))
            if len(idx)>n_samp:
                idx = np.sort(np.random.choice(idx, size=n_samp, replace=False))
            mask[idx] = True
        bin_str = 'bin_' + str(bin_index)
        mask_dict[bin_str] = mask.copy()
        mask_all |= mask

    cat['pz_bin'] = 0
    for bin_index in range(1, 5):
        cat['pz_bin'][mask_dict['bin_{}'.format(bin_index)]] = bin_index

    # # Sanity checks: there should be no overlap between bins
    # print(np.sum(mask_dict['bin_1'] & mask_dict['bin_2']))
    # print(np.sum(mask_dict['bin_1'] & mask_dict['bin_4']))
    # print(np.sum(mask_dict['bin_3'] & mask_dict['bin_4']))
    # print(np.sum(mask_dict['bin_2'] & mask_dict['bin_3']))
    # print(np.sum(mask_dict['bin_1'] | mask_dict['bin_2'] | mask_dict['bin_3'] | mask_dict['bin_4']))
    # print(np.sum(mask_dict['bin_1'])+np.sum(mask_dict['bin_2'])+np.sum(mask_dict['bin_3'])+np.sum(mask_dict['bin_4']))

    # plt.figure(figsize=(15, 10))
    # for bin_index in range(1, 5):
    #     plt.hist(cat['Z'][mask_dict['bin_{}'.format(bin_index)]], bins=tt['zmin'], weights=1/y1_area*np.ones(np.sum(mask_dict['bin_{}'.format(bin_index)])), alpha=0.5, color='C'+str(bin_index))
    #     plt.plot((tt['zmin']+tt['zmax'])/2, tt['bin_{}_combined'.format(bin_index)], color='C'+str(bin_index))
    # plt.grid(alpha=0.5)
    # plt.show()


    cat_all = cat.copy()

    n_processes = 1
    npix = hp.nside2npix(nside)


    def get_z_variation(pix_idx):

        pix_list = pix_unique[pix_idx]

        hp_table = Table()
        hp_table['HPXPIXEL'] = pix_list
        hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_list, nest=False, lonlat=True)
        bin_str = 'bin_'+str(pz_bin)+'_'
        for col in [bin_str+'z_mean', bin_str+'z_median', bin_str+'z_l68', bin_str+'z_u68', bin_str+'z_l95', bin_str+'z_u95', bin_str+'z_nmad']:
            hp_table[col] = np.zeros(len(hp_table), dtype=float)
        hp_table[bin_str+'n_objects'] = np.zeros(len(hp_table), dtype=int)
        for index in np.arange(len(pix_idx)):
            idx = pixorder[pixcnts[pix_idx[index]]:pixcnts[pix_idx[index]+1]]
            hp_table[bin_str+'z_mean'][index] = np.mean(cat['Z'][idx])
            hp_table[bin_str+'z_median'][index], hp_table[bin_str+'z_l68'][index], hp_table[bin_str+'z_u68'][index], hp_table[bin_str+'z_l95'][index], hp_table[bin_str+'z_u95'][index] =\
                np.percentile(cat['Z'][idx], [50, 16., 84., 2.5, 97.5])
            hp_table[bin_str+'z_nmad'][index] = nmad(cat['Z'][idx])
            hp_table[bin_str+'n_objects'][index] = len(idx)

        return hp_table


    # z variation
    for index, pz_bin in enumerate(range(1, 5)):
        print(pz_bin)
        mask = (cat_all['pz_bin']==pz_bin)
        cat = cat_all[mask].copy()

        pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'], cat['DEC'], lonlat=True)
        pix_unique, pixcnts = np.unique(pix_allobj, return_counts=True)

        pixcnts = np.insert(pixcnts, 0, 0)
        pixcnts = np.cumsum(pixcnts)

        pixorder = np.argsort(pix_allobj)

        pix_idx_split = np.array_split(np.arange(len(pix_unique)), n_processes)

        # start multiple worker processes
        with Pool(processes=n_processes) as pool:
            res = pool.map(get_z_variation, pix_idx_split)

        if index==0:
            hp_table = vstack(res)
        else:
            tmp = vstack(res)
            tmp.remove_columns(['RA', 'DEC'])
            hp_table = join(hp_table, tmp, keys='HPXPIXEL', join_type='outer')

    hp_table = hp_table.filled(0)
    hp_table.sort('HPXPIXEL')
    hp_table.write('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/misc/lrg_pzbins_specz_stats_{}-AbacusSummit_Y1v1_mock{}.fits'.format(nside, mock_subset), overwrite=True)

