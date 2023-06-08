# The decam2014 filter curves were used in iron reduction
# Compute no-reddening and (SFD-)reddened magnitudes in BASS/MzLS and DECam filters
# iron tag: "source /global/common/software/desi/desi_environment.sh 23.1"
# desispec/0.56.5
# speclite/0.16

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from multiprocessing import Pool

import speclite.filters
from desispec.magnitude import compute_ab_mag
from desiutil.dust import dust_transmission

fluxunits = 1e-17  # units.erg / units.s / units.cm**2 / units.Angstrom


def load_legacy_survey_filter(band, photsys):
    """
    Uses speclite.filters to load the filter transmission
    Returns speclite.filters.FilterResponse object

    Args:
        band: filter pass-band in "G","R","Z","W1","W2"
        photsys: "N" or "S" for North (BASS+MzLS) or South (CTIO/DECam)
    """
    filternamemap=None
    if band[0].upper()=="W":  # it's WISE
        filternamemap = "wise2010-{}".format(band.upper())
    elif band.upper() in ['G', 'R', 'I', 'Z']:
        if photsys=="N":
            if band.upper() in ["G", "R"]:
                filternamemap="BASS-{}".format(band.lower())
            else:
                filternamemap="MzLS-z"
        elif photsys=="S":
            ##################################################
            filternamemap="decam2014-{}".format(band.lower())
            ##################################################
        else:
            raise ValueError("unknown photsys '{}', known ones are 'N' and 'S'".format(photsys))
    else:
        raise ValueError("unknown band '{}', known ones are 'G','R', 'I', 'Z','W1' and 'W2'".format(photsys))

    filter_response = speclite.filters.load_filter(filternamemap)
    return filter_response


filter_curves = dict()
for band in ['G', 'R', 'I', 'Z']:
    for photsys in ['N', 'S']:
        filtername = band + photsys
        if filtername=='IN':  # No i band in North
            continue
        filter_curves[filtername] = load_legacy_survey_filter(band=band, photsys=photsys)
for band in ['W1', 'W2']:
    filter_curves[band] = load_legacy_survey_filter(band=band, photsys=None)

exp = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/iron/exposures-iron.fits', ext='EXPOSURES'))
print(len(exp))


def get_std(exp_index):

    cat_stack = []

    for petal_loc in range(10):

        night = exp['NIGHT'][exp_index]
        expid = exp['EXPID'][exp_index]
        expid_str = str(expid).zfill(8)
        fn = f'/global/cfs/cdirs/desi/spectro/redux/iron/exposures/{night}/{expid_str}/stdstars-{petal_loc}-{expid_str}.fits.gz'

        if not os.path.isfile(fn):
            continue

        cat1 = Table(fitsio.read(fn, ext='METADATA'))
        cat2 = Table(fitsio.read(fn, ext='FIBERMAP'))
        # cat1.sort('TARGETID')
        # cat2.sort('TARGETID')
        assert np.all(cat1['TARGETID']==cat2['TARGETID']) and np.all(cat1['FIBER']==cat2['FIBER'])

        cat2.remove_columns(['TARGETID', 'FIBER'])
        cat = hstack([cat1, cat2])
        # print(expid, len(cat))

        stdwave = fitsio.read(fn, ext='WAVELENGTH')
        # print(stdwave.shape)
        models = fitsio.read(fn, ext='FLUX')
        # print(models.shape)

        for band in ['G', 'R', 'I', 'Z']:
            for photsys in ['N', 'S']:
                if band+photsys=='IN':  # No i band in North
                    continue
                cat['MODEL_{}MAG_{}_REDDENED'.format(band, photsys)] = 0.
                cat['MODEL_{}MAG_{}'.format(band, photsys)] = 0.
                cat['MODEL_{}MAG_{}_REDDENED_SPECLITE'.format(band, photsys)] = 0.
                cat['MODEL_{}MAG_{}_SPECLITE'.format(band, photsys)] = 0.
        for band in ['W1', 'W2']:
            cat['MODEL_{}MAG_REDDENED'.format(band)] = 0.
            cat['MODEL_{}MAG'.format(band)] = 0.
            cat['MODEL_{}MAG_REDDENED_SPECLITE'.format(band)] = 0.
            cat['MODEL_{}MAG_SPECLITE'.format(band)] = 0.

        for index in range(len(cat)):

            ebv = cat['EBV'][index]
            # photsys = cat['PHOTSYS'][index]

            model = models[index].copy()
            model_no_reddening = model / dust_transmission(stdwave, ebv)  # undo the dereddening

            for band in ['G', 'R', 'I', 'Z']:
                for photsys in ['N', 'S']:
                    filtername = band + photsys
                    if filtername=='IN':  # No i band in North
                        continue
                    filt_ww, filt_tt = filter_curves[filtername].wavelength.copy(), filter_curves[filtername].response.copy()
                    # Reddened model magnitude using Julien's script
                    cat['MODEL_{}MAG_{}_REDDENED'.format(band, photsys)][index] = compute_ab_mag(stdwave, model, filt_ww, filt_tt)
                    # Zero-reddening model magnitude using Julien's script
                    cat['MODEL_{}MAG_{}'.format(band, photsys)][index] = compute_ab_mag(stdwave, model_no_reddening, filt_ww, filt_tt)
                    # Reddened model magnitude using speclite (used in iron reduction)
                    cat['MODEL_{}MAG_{}_REDDENED_SPECLITE'.format(band, photsys)][index] = filter_curves[filtername].get_ab_magnitude(model * fluxunits, stdwave.copy())
                    # Zero-reddening model magnitude using speclite (used in iron reduction)
                    cat['MODEL_{}MAG_{}_SPECLITE'.format(band, photsys)][index] = filter_curves[filtername].get_ab_magnitude(model_no_reddening * fluxunits, stdwave.copy())
            for band in ['W1', 'W2']:
                filt_ww, filt_tt = filter_curves[band].wavelength.copy(), filter_curves[band].response.copy()
                # Reddened model magnitude using Julien's script
                cat['MODEL_{}MAG_REDDENED'.format(band)][index] = compute_ab_mag(stdwave, model, filt_ww, filt_tt)
                # Zero-reddening model magnitude using Julien's script
                cat['MODEL_{}MAG'.format(band)][index] = compute_ab_mag(stdwave, model_no_reddening, filt_ww, filt_tt)
                # Reddened model magnitude using speclite (used in iron reduction)
                cat['MODEL_{}MAG_REDDENED_SPECLITE'.format(band)][index] = filter_curves[band].get_ab_magnitude(model * fluxunits, stdwave.copy())
                # Zero-reddening model magnitude using speclite (used in iron reduction)
                cat['MODEL_{}MAG_SPECLITE'.format(band)][index] = filter_curves[band].get_ab_magnitude(model_no_reddening * fluxunits, stdwave.copy())

        cat_stack.append(cat)

    if len(cat_stack)==0:
        return None
    else:
        cat_stack = vstack(cat_stack)
        cat_stack['EXPID'] = expid
        return cat_stack


n_processess = 128
with Pool(processes=n_processess) as pool:
    res = pool.map(get_std, np.arange(len(exp)))

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

cat = vstack(res)
cat.write('/pscratch/sd/r/rongpu/ebv/desi_std/desi_standard_stars_iron_add_magnitudes_decam2014.fits')

