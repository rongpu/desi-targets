from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack
import fitsio

from multiprocessing import Pool
import healpy as hp


# resolve = 'noresolve'
resolve = 'resolve'

field = str(sys.argv[1])
field = field.lower()

if field=='south':
    photsys = 'S'
elif field=='north':
    photsys = 'N'

min_nobs = 1
# maskbits = sorted([1, 13])
# maskbits = sorted([1, 8, 9, 11, 12, 13])
maskbits = sorted([1, 11, 12, 13])

n_processes = 32

n_randoms_catalogs = 32  # There are 200 random catalogs in total

nsides = [64, 128, 256]
randoms_columns = ['RA', 'DEC', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'MASKBITS', 'PHOTSYS']

if resolve=='resolve':
    randoms_paths = sorted(glob.glob('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-[0-9]*.fits'))
elif resolve=='noresolve':
    randoms_paths = sorted(glob.glob('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/noresolve/{}/randoms-noresolve-*.fits'.format(field)))
print(len(randoms_paths))

randoms_paths = randoms_paths[:n_randoms_catalogs]
print(len(randoms_paths))

randoms_density = 2500

output_dir = '/global/cfs/cdirs/desi/users/rongpu/imaging_sys/randoms_stats/0.49.0/{}/counts'.format(resolve)


def apply_mask(randoms, min_nobs, maskbits):

    mask = (randoms['NOBS_G']>=min_nobs) & (randoms['NOBS_R']>=min_nobs) & (randoms['NOBS_Z']>=min_nobs)

    mask_clean = np.ones(len(randoms), dtype=bool)
    for bit in maskbits:
        mask_clean &= (randoms['MASKBITS'] & 2**bit)==0
    # print(np.sum(~mask_clean)/len(mask_clean))

    mask &= mask_clean

    return mask


def count_randoms(randoms_path):

    print(randoms_path)
    if resolve=='resolve':
        randoms_index_str = os.path.basename(randoms_path).replace('randoms-', '').replace('.fits', '')
    elif resolve=='noresolve':
        randoms_index_str = os.path.basename(randoms_path).replace('randoms-noresolve-', '').replace('.fits', '')

    randoms = Table(fitsio.read(randoms_path, columns=randoms_columns))
    # print(len(randoms))

    if fitsio.read_header(randoms_path, ext=1)['DENSITY']!=randoms_density:
        raise ValueError

    if resolve=='resolve':
        mask = randoms['PHOTSYS']==photsys
        randoms = randoms[mask]

    mask = apply_mask(randoms, min_nobs, maskbits)
    randoms = randoms[mask]

    for nside in nsides:
        npix = hp.nside2npix(nside)

        # pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)
        # print('Healpix size = {:.5f} sq deg'.format(pix_area))

        pix = hp.pixelfunc.ang2pix(nside, randoms['RA'], randoms['DEC'], lonlat=True)
        pix_unique, pix_count = np.unique(pix, return_counts=True)
        pix_count_all = np.zeros(npix, dtype=int)
        pix_count_all[pix_unique] = pix_count

        output_path = os.path.join(output_dir, 'minobs_{}_maskbits_{}'.format(min_nobs, ''.join([str(tmp) for tmp in maskbits])), '{}_nside_{}_minobs_{}_maskbits_{}_{}.npy'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits]), randoms_index_str))
        if not os.path.isdir(os.path.dirname(output_path)):
            try:
                os.makedirs(os.path.dirname(output_path))
            except:
                pass

        np.save(output_path, pix_count_all)

    return None


if __name__ == '__main__':

    print('Start!')

    time_start = time.time()

    # # start multiple worker processes
    # with Pool(processes=n_processes) as pool:
    #     pool.map(count_randoms, randoms_paths)
   
    for randoms_path in randoms_paths:
        count_randoms(randoms_path)

    print(time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))

    # Combine the results into a single table

    for nside in nsides:

        print(nside)

        npix = hp.nside2npix(nside)
        pix_area = hp.pixelfunc.nside2pixarea(nside, degrees=True)

        hp_table = Table()
        hp_table['hp_idx'] = np.arange(npix)
        hp_table['ra'], hp_table['dec'] = hp.pixelfunc.pix2ang(nside, hp_table['hp_idx'], nest=False, lonlat=True)
        hp_table['count'] = 0

        output_paths = sorted(glob.glob(os.path.join(output_dir, 'minobs_{}_maskbits_{}'.format(min_nobs, ''.join([str(tmp) for tmp in maskbits])), '{}_nside_{}_minobs_{}_maskbits_{}_*.npy'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])))))
        print(len(output_paths))

        for output_path in output_paths:
            count = np.load(output_path)
            hp_table['count'] += count

        total_randoms_density = randoms_density * len(output_paths)
        hp_table['pix_frac'] = hp_table['count']/(total_randoms_density*pix_area)

        hp_table.write((os.path.join(output_dir, 'counts_{}_nside_{}_minobs_{}_maskbits_{}.fits'.format(field, nside, min_nobs, ''.join([str(tmp) for tmp in maskbits])))))
