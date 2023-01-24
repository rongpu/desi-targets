# srun -N 1 python rf_photoz_with_magnification.py south 0.99

from __future__ import division, print_function
import sys, os, warnings, gc, time, glob
from pathlib import Path
from multiprocessing import Pool

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, hstack, vstack
import fitsio

import joblib

from desitarget.targets import decode_targetid, encode_targetid

###############
test_run = False
###############

n_process = 128

n_folds = 10

n_estimators = 100  # Number of trees in a forest
n_perturb = 20  # Number of perturbed sample

pz_dtype = [('Z_PHOT_MEAN', 'f4'), ('Z_PHOT_MEDIAN', 'f4'), ('Z_PHOT_STD', 'f4'),
            ('Z_PHOT_L68', 'f4'), ('Z_PHOT_U68', 'f4'), ('Z_PHOT_L95', 'f4'), ('Z_PHOT_U95', 'f4')]

columns = ['RELEASE', 'BRICKID', 'OBJID', 'TYPE', 'RA', 'DEC', 'DCHISQ', 'EBV',
'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_R',
'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2',
'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z', 'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2',
'NOBS_G', 'NOBS_R', 'NOBS_Z',
'SHAPE_R', 'SHAPE_R_IVAR', 'SHAPE_E1', 'SHAPE_E2']

####################################################################################################################


def unpack(cat):
    gmag = 22.5-2.5*np.log10(cat['FLUX_G_EC'])
    rmag = 22.5-2.5*np.log10(cat['FLUX_R_EC'])
    zmag = 22.5-2.5*np.log10(cat['FLUX_Z_EC'])
    w1mag = 22.5-2.5*np.log10(cat['FLUX_W1_EC'])
    w2mag = 22.5-2.5*np.log10(cat['FLUX_W2_EC'])
    rfibermag = 22.5-2.5*np.log10(cat['FIBERFLUX_R_EC'])
    radius1 = np.array(cat['SHAPE_R'])
    return gmag, rmag, zmag, w1mag, w2mag, rfibermag, radius1


def compute_photoz(split_idx, kf_index, magnification):

    if test_run:
        split_idx = split_idx.copy()
        split_idx = split_idx[::100]
        print(len(split_idx))

    if len(split_idx)==0:
        return None

    cat = cat_all[split_idx].copy()

    cat_id_full = cat[['RELEASE', 'BRICKID', 'OBJID']].copy()

    # Sanity check to make sure that there are not extra spaces in type names
    types = ['PSF', 'REX', 'EXP', 'DEV', 'SER', 'DUP']
    for tmp in np.unique(cat['TYPE']):
        if tmp not in types:
            raise ValueError('Type {} not recognized'.format(tmp))

    full_size = len(cat)

    mask_bad = np.full(len(cat), False, dtype=bool)
    mask_bad |= ~((cat['NOBS_G']>=1) & (cat['NOBS_R']>=1) & (cat['NOBS_Z']>=1))
    mask_bad |= ~((cat['FLUX_IVAR_G']>0) & (cat['FLUX_IVAR_R']>0) & (cat['FLUX_IVAR_Z']>0))
    mask_bad |= (cat['TYPE']=='DUP')

    print('{:} ({:.1f}%) objects removed'.format(np.sum(mask_bad), np.sum(mask_bad)/len(mask_bad)*100))

    if np.sum(~mask_bad)==0:
        return None

    cat = cat[~mask_bad]
    print(len(cat))

    # Assign inf to invalid FLUX_IVAR_W1 and FLUX_IVAR_W2
    # (Not sure if it actually happens for TYPE!='DUP' objects though)
    mask = cat['FLUX_IVAR_W1']==0
    cat['FLUX_IVAR_W1'][mask] = np.inf
    mask = cat['FLUX_IVAR_W2']==0
    cat['FLUX_IVAR_W2'][mask] = np.inf

    # Assign inf to invalid zero SHAPE_R_IVAR
    mask = cat['SHAPE_R_IVAR']==0
    cat['SHAPE_R_IVAR'][mask] = np.inf

    # axis ratio
    e = np.array(np.sqrt(cat['SHAPE_E1']**2+cat['SHAPE_E2']**2))
    q = (1+e)/(1-e)

    # shape probability (definition of shape probability in Soo et al. 2017)
    p = np.ones(len(cat))*0.5
    # DCHISQ[:, 2] is DCHISQ_EXP; DCHISQ[:, 3] is DCHISQ_DEV
    mask_chisq = (cat['DCHISQ'][:, 3]>0) & (cat['DCHISQ'][:, 2]>0)
    p[mask_chisq] = cat['DCHISQ'][:, 3][mask_chisq]/(cat['DCHISQ'][:, 3]+cat['DCHISQ'][:, 2])[mask_chisq]

    ####################################################################################################################

    cat_pz = Table(data=np.full(len(cat), -99, dtype=pz_dtype))

    # print('Computing photo-z\'s and photo-z errors')

    col_list = ['FLUX_G_EC', 'FLUX_R_EC', 'FLUX_Z_EC', 'FLUX_W1_EC', 'FLUX_W2_EC', 'FIBERFLUX_R_EC']
    mag_max = 30
    mag_fill = 100

    dtype1 = [('FLUX_G_EC', 'f4'), ('FLUX_R_EC', 'f4'), ('FLUX_Z_EC', 'f4'), ('FLUX_W1_EC', 'f4'), ('FLUX_W2_EC', 'f4'), ('SHAPE_R', 'f4')]
    cat1 = Table(data=np.zeros(len(cat), dtype=dtype1))

    z_phot_array = np.zeros((len(cat1), n_perturb*n_estimators), dtype='float32')

    for tree_index in range(n_estimators):

        # print(tree_index*n_perturb, '/', n_estimators*n_perturb)

        # Predict!
        for perturb_index in range(n_perturb):

            cat1['FLUX_G_EC'] = magnification * np.array((cat['FLUX_G']+np.random.randn(len(cat))/np.sqrt(cat['FLUX_IVAR_G']))/cat['MW_TRANSMISSION_G'], dtype='float32')
            # cat1['FLUX_R_EC'] = magnification * np.array((cat['FLUX_R']+np.random.randn(len(cat))/np.sqrt(cat['FLUX_IVAR_R']))/cat['MW_TRANSMISSION_R'], dtype='float32')
            cat1['FLUX_Z_EC'] = magnification * np.array((cat['FLUX_Z']+np.random.randn(len(cat))/np.sqrt(cat['FLUX_IVAR_Z']))/cat['MW_TRANSMISSION_Z'], dtype='float32')
            cat1['FLUX_W1_EC'] = magnification * np.array((cat['FLUX_W1']+np.random.randn(len(cat))/np.sqrt(cat['FLUX_IVAR_W1']))/cat['MW_TRANSMISSION_W1'], dtype='float32')
            cat1['FLUX_W2_EC'] = magnification * np.array((cat['FLUX_W2']+np.random.randn(len(cat))/np.sqrt(cat['FLUX_IVAR_W2']))/cat['MW_TRANSMISSION_W2'], dtype='float32')
            cat1['SHAPE_R'] = magnification * np.array(cat['SHAPE_R']+np.random.randn(len(cat))/np.sqrt(cat['SHAPE_R_IVAR']), dtype='float32')

            noise_flux_r = np.random.randn(len(cat))
            cat1['FLUX_R_EC'] = magnification * np.array((cat['FLUX_R']+noise_flux_r/np.sqrt(cat['FLUX_IVAR_R']))/cat['MW_TRANSMISSION_R'], dtype='float32')
            cat1['FIBERFLUX_R_EC'] = magnification * np.array((cat['FIBERFLUX_R']+cat['FIBERFLUX_R']/cat['FLUX_R']*noise_flux_r/np.sqrt(cat['FLUX_IVAR_R']))/cat['MW_TRANSMISSION_R'], dtype='float32')

            # Fill in negative fluxes
            for index in range(len(col_list)):
                mask = (cat1[col_list[index]]<10**(0.4*(22.5-mag_max))) | (~np.isfinite(cat1[col_list[index]]))
                cat1[col_list[index]][mask] = 10**(0.4*(22.5-mag_fill))

            gmag1, rmag1, zmag1, w1mag1, w2mag1, rfibermag1, radius1 = unpack(cat1)
            data1 = np.column_stack((gmag1-rmag1, rmag1-zmag1, zmag1-w1mag1, w1mag1-w2mag1, rmag1, rfibermag1, radius1, q, p))

            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=UserWarning)
                z_phot_array[:, tree_index*n_perturb+perturb_index] = regrf_all[tree_index].predict(data1)

        # clear cache
        gc.collect()

    cat_pz['Z_PHOT_MEAN'] = np.mean(z_phot_array, axis=1)
    cat_pz['Z_PHOT_STD'] = np.std(z_phot_array, axis=1)
    # WARNING: z_phot_array becomes undefined after this
    cat_pz['Z_PHOT_L95'], cat_pz['Z_PHOT_L68'], cat_pz['Z_PHOT_MEDIAN'], cat_pz['Z_PHOT_U68'], cat_pz['Z_PHOT_U95'] = \
        np.percentile(z_phot_array, [2.5, 16., 50., 84., 97.5], axis=1, overwrite_input=True)

    cat_pz_full = Table(data=np.full(full_size, -99, dtype=pz_dtype))
    cat_pz_full[~mask_bad] = cat_pz

    # Add RELEASE, BRICKID and OBJID
    cat_pz_full = hstack([cat_id_full, cat_pz_full])

    cat_pz_full['KFOLD'] = np.full(len(cat_pz_full), kf_index, dtype=np.int16)
    mask = cat_pz_full['Z_PHOT_MEDIAN']==-99
    cat_pz_full['KFOLD'][mask] = -99

    return cat_pz_full


if __name__ == '__main__':

    time_start = time.time()

    field = str(sys.argv[1])
    magnification = float(sys.argv[2])
    print(field, magnification)

    if (field!='south') and (field!='north'):
        raise ValueError('field must be either \"north\" or \"south\"')

    output_path = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/extended_lrg_magnification_pz_{}_{:g}.fits'.format(field, magnification)
    chunks_dir = '/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/tmp/kf_chunks/'+field

    cat_all = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/lrg_xcorr/magnification/extended_lrg_magnification_{}.fits'.format(field), columns=columns))
    cat_size = len(cat_all)

    np.random.seed(1456+cat_size)
    kf_idx = np.random.choice(n_folds, size=cat_size)

    cat_all['TARGETID'] = encode_targetid(cat_all['OBJID'], cat_all['BRICKID'], cat_all['RELEASE'])

    # Full spec-z truth catalog (not all objects were used for training):
    specz_full_path = '/global/cfs/cdirs/desi/users/rongpu/data/ls_dr9.0_desi_photoz/truth/truth_combined_dr9.0_desi_20221021_{}.fits'.format(field)
    # Catalog used for training
    specz_train_path = '/global/cfs/cdirs/desi/users/rongpu/data/ls_dr9.0_desi_photoz/truth/truth_combined_dr9.0_desi_20221021_{}_training.fits'.format(field)

    specz_train = Table(fitsio.read(specz_train_path, columns=['TARGETID', 'redshift', 'survey', 'kfold']))
    # specz_train['survey'] = np.array(specz_train['survey']).astype('str')
    train_tid_all = np.array(specz_train['TARGETID'])

    # specz_full = Table(fitsio.read(specz_full_path, columns=['TARGETID', 'redshift', 'survey']))
    # # specz_full['survey'] = np.array(specz_full['survey']).astype('str')

    print('Start!')

    for kf_index in range(n_folds):

        chunk_output_path = os.path.join(chunks_dir, 'extended_lrg_magnification_pz_{}_{:g}.fits'.format(kf_index+1, magnification))
        if os.path.isfile(chunk_output_path):
            continue
        Path(chunk_output_path).touch(exist_ok=False)

        print('K-fold', kf_index+1)

        mask = specz_train['kfold']==kf_index
        test_tid = np.array(specz_train['TARGETID'][mask])

        mask_kf = (kf_idx==kf_index) & (~np.in1d(cat_all['TARGETID'], train_tid_all))  # randomly assigned KF membership
        mask_train = np.in1d(cat_all['TARGETID'], test_tid)  # objects in the spec-z KF test set
        cat_idx = np.where(mask_kf | mask_train)[0]

        # split among the processes
        chunk_size = float(2e4)
        n_split = int(np.ceil(len(cat_idx)/chunk_size))
        split_idx_list = np.array_split(cat_idx, n_split)
        if n_folds==0:
            print('n_split', n_split)

        tree_dir = '/global/cfs/cdirs/desi/users/rongpu/data/ls_dr9.0_desi_photoz/individual_trees/20221021_10-fold-{}/fold-{}'.format(field, kf_index+1)

        # Load all the single pre-trained trees
        regrf_all = []
        for tree_index in range(n_estimators):
            regrf_all.append(joblib.load(os.path.join(tree_dir, 'regrf_20221021_{:d}.pkl'.format(tree_index))))

        zipped_arg_list = list(zip(split_idx_list, [kf_index]*n_split, [magnification]*n_split))

        # start multiple worker processes
        with Pool(processes=n_process) as pool:
            res = pool.starmap(compute_photoz, zipped_arg_list)

        # Remove None elements from the list
        for index in range(len(res)-1, -1, -1):
            if res[index] is None:
                res.pop(index)

        cat = vstack(res)
        cat.write(chunk_output_path, overwrite=True)

        print('K-fold {} done!'.format(kf_index+1), time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))

    print('All done!', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))

    pz_stack = []
    for kf_index in range(n_folds):
        chunk_path = os.path.join(chunks_dir, 'extended_lrg_magnification_pz_{}_{:g}.fits'.format(kf_index+1, magnification))
        if os.stat(chunk_path).st_size!=0:
            pz_stack.append(Table(fitsio.read(chunk_path)))
    pz = vstack(pz_stack)

    pz['TARGETID'] = encode_targetid(pz['OBJID'], pz['BRICKID'], pz['RELEASE'])

    # Here matching tt to cat_all
    if len(cat_all)!=len(pz) or not np.all(np.unique(cat_all['TARGETID'])==np.unique(pz['TARGETID'])):
        raise ValueError
    t1_reverse_sort = np.array(cat_all['TARGETID']).argsort().argsort()
    pz = pz[np.argsort(pz['TARGETID'])[t1_reverse_sort]]

    # Double check
    if len(cat_all)!=len(pz) or not np.all(np.unique(cat_all['TARGETID'])==np.unique(pz['TARGETID'])):
        raise ValueError

    pz.remove_column('TARGETID')

    pz.write(output_path)
