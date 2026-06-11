#!/usr/bin/envpython

# To run:
# source /dvs_ro/cfs/cdirs/desi/software/desi_environment.sh master
# salloc -N 1 -C haswell -t 01:00:00 --qos interactive -L SCRATCH,project
# python desi-per-tracer-allspec.py --tracer LRG --numproc 32
# or for SV3: python desi-per-tracer-allspec.py --survey sv3 --tracer LRG --numproc 32

from astropy.io import fits
from astropy.table import Table, vstack
import fitsio
import numpy as np
import os
from glob import glob
import sys
from desitarget.geomask import match
from desitarget.internal import sharedmem
# from raichoorlib import get_foii,get_date
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument(
    "--survey",
    help="survey (default=main)",
    type=str,
    default="main",
    choices=["sv3", "main"]
)
parser.add_argument(
    "--tracer",
    help="tracer",
    type=str,
    default=None,
    choices=["MWS_ANY", "BGS_ANY", "LRG", "ELG", "QSO"]
)
parser.add_argument(
    "--specprod",
    help="daily (default=daily)",
    type=str,
    default="daily",
)
parser.add_argument(
    "--numproc",
    help="number of concurrent processes to use; set to 0 to not process (default=0)",
    type=int,
    default=1,
)
args = parser.parse_args()
for kwargs in args._get_kwargs():
    print(kwargs)


#outdir = "{}/{}/".format(os.getenv("CSCRATCH"), args.survey)
# HACK ===========
outdir = os.path.join('/global/cfs/cdirs/desi/users/rongpu/data/', args.survey)
# HACK ===========

pipedir = "{}/spectro/redux/{}/tiles/cumulative".format(os.getenv("DESI_ROOT"), args.specprod)
tsp_separator = ","

if args.survey == "sv3":
    from desitarget.sv3.sv3_targetmask import desi_mask
    dtkey = "SV3_DESI_TARGET"
    elg_tracers = ["ELG", "ELG_HIP", "ELG_LOP"]
if args.survey == "main":
    from desitarget.targetmask import desi_mask
    dtkey = "DESI_TARGET"
    elg_tracers = ["ELG", "ELG_LOP", "ELG_HIP", "ELG_VLO"]

# columns we keep
keys = {}
keys["ZBEST"] = ["TARGETID", "CHI2", "COEFF", "Z", "ZERR", "ZWARN", "NPIXELS", "SPECTYPE", "SUBTYPE", "NCOEFF", "DELTACHI2", "NUMEXP", "NUMTILE"]
keys["FIBERMAP"] = ["PETAL_LOC", "DEVICE_LOC", "LOCATION", "FIBER", "FIBERSTATUS", "TARGET_RA", "TARGET_DEC", "PMRA", "PMDEC", "REF_EPOCH", "LAMBDA_REF", "FA_TARGET", "FA_TYPE", "OBJTYPE", "FIBERASSIGN_X", "FIBERASSIGN_Y", "PRIORITY", "SUBPRIORITY", "OBSCONDITIONS", "RELEASE", "BRICKID", "BRICK_OBJID", "MORPHTYPE", "FLUX_G", "FLUX_R", "FLUX_Z", "FLUX_IVAR_G", "FLUX_IVAR_R", "FLUX_IVAR_Z", "MASKBITS", "REF_ID", "REF_CAT", "GAIA_PHOT_G_MEAN_MAG", "GAIA_PHOT_BP_MEAN_MAG", "GAIA_PHOT_RP_MEAN_MAG", "PARALLAX", "BRICKNAME", "EBV", "FLUX_W1", "FLUX_W2", "FLUX_IVAR_W1", "FLUX_IVAR_W2", "FIBERFLUX_G", "FIBERFLUX_R", "FIBERFLUX_Z", "FIBERTOTFLUX_G", "FIBERTOTFLUX_R", "FIBERTOTFLUX_Z", "SERSIC", "SHAPE_R", "SHAPE_E1", "SHAPE_E2", "PHOTSYS", "PRIORITY_INIT", "NUMOBS_INIT", "TILEID"]
if args.survey == "sv3":
    keys["FIBERMAP"] += ["SV3_DESI_TARGET", "SV3_BGS_TARGET", "SV3_MWS_TARGET", "SV3_SCND_TARGET"]
if args.survey == "main":
    keys["FIBERMAP"] += ["DESI_TARGET", "BGS_TARGET", "MWS_TARGET", "SCND_TARGET"]
# own + scores
okeys = ["FN", "HASCOADD", "THRUNIGHT",
    # "FOII", "FOII_ERR",
    # "GTOT", "GFIB", "GR", "RZ", "RW1", "RW2",
    "TSNR2_ELG", "TSNR2_LRG", 'TSNR2_BGS',
    "COADD_EFFTIME_SPEC", "COADD_EFFTIME_ETC", "COADD_EFFTIME_GFA",
    "ZDONE"]
ofmts = ["A200", "L", "K", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "E", "L", "E", "E", "E"]


# tsnr
tsnr = fits.open("{}/spectro/redux/{}/tsnr-exposures.fits".format(os.getenv("DESI_ROOT"), args.specprod))[1].data

# tiles
status_tiles = Table.read("{}/spectro/redux/{}/tiles.csv".format(os.getenv("DESI_ROOT"), args.specprod))

# tiles
tiles = Table.read("{}/ops/tiles-{}.ecsv".format(os.getenv("DESI_SURVEYOPS"), args.survey))
if args.tracer[:3] in ["LRG", "ELG", "QSO"]:
    tiles = tiles[tiles["PROGRAM"] == "DARK"]
if args.tracer[:3] in ["BGS", "MWS"]:
    tiles = tiles[tiles["PROGRAM"] == "BRIGHT"]
tileids = tiles["TILEID"]
# grabbing the latest night for each tileid
# and restricting to observed tileids
thrunights = np.array(["00000000" for i in tileids])
for i in range(len(tiles)):
    # require 20* to reject "stash"
    fns = np.sort(glob("{}/{}/20*".format(pipedir, tiles["TILEID"][i])))
    if len(fns) > 0:
        thrunights[i] = fns[-1].split("/")[-1]
keep = thrunights != "00000000"
tileids, thrunights = tileids[keep], thrunights[keep]


# output fits
outfits = "{}/{}-{}-{}-thru{}.fits".format(outdir, args.survey, args.tracer.lower(), args.specprod, np.unique(thrunights)[-1])





#
def encode_tsp(tileid, thrunight, petal):
    return tsp_separator.join(["{}".format(tileid), thrunight, "{}".format(petal)])

#
def decode_tsp(tileid_thrunight_petal):
    tileid = int(tileid_thrunight_petal.split(tsp_separator)[0])
    thrunight = tileid_thrunight_petal.split(tsp_separator)[1]
    petal = int(tileid_thrunight_petal.split(tsp_separator)[2])
    return tileid, thrunight, petal


# zbest files
def get_zbfn(tileid, thrunight, petal):
    return "{}/{}/{}/zbest-{}-{}-thru{}.fits".format(
                pipedir, tileid, thrunight, petal, tileid, thrunight)

# 
def get_outfn(tileid, thrunight, petal):
    return "{}/{}-per-petal/{}-{}-{}-{}-thru{}-{}.fits".format(
                        outdir, args.survey, args.survey, args.tracer.lower(), args.specprod, tileid, thrunight, petal)

#
# def get_magcol(d):
#     exts = {"g":3.214, "r":2.165, "z":1.211, "w1":0.184, "w2":0.113}
#     gtot = 22.5-2.5*np.log10(d["flux_g"]) - exts["g"]*d["ebv"]
#     gfib = 22.5-2.5*np.log10(d["fiberflux_g"]) - exts["g"]*d["ebv"]
#     gr   = -2.5*np.log10(d["flux_g"]/d["flux_r"]) - (exts["g"]-exts["r"])*d["ebv"]
#     rz   = -2.5*np.log10(d["flux_r"]/d["flux_z"]) - (exts["r"]-exts["z"])*d["ebv"]
#     rw1  = -2.5*np.log10(d["flux_r"]/d["flux_w1"]) - (exts["r"]-exts["w1"])*d["ebv"]
#     rw2  = -2.5*np.log10(d["flux_r"]/d["flux_w2"]) - (exts["r"]-exts["w2"])*d["ebv"]
#     return gtot, gfib, gr, rz, rw1, rw2

# efftime per reduction
def get_efftime(zbfn):
    expids = np.unique(fits.open(zbfn)["FIBERMAP"].data["EXPID"])
    keep = np.in1d(tsnr["EXPID"], expids)
    efftime_spec = tsnr["EFFTIME_SPEC"][keep].sum()
    efftime_etc = tsnr["EFFTIME_ETC"][keep].sum()
    efftime_gfa = tsnr["EFFTIME_GFA"][keep].sum()
    return efftime_spec, efftime_etc, efftime_gfa


# weights to account for elgxqso over-representation
# per-coadd
# (n_elgqso * w) / (n_elgqso * w + n_elg-n_elgqso) = fracp
# n_elgqso * w = fracp * (n_elgqso * w + n_elg-n_elgqso)
# (n_elgqso - fracp * n_elgqso) * w = fracp * (n_elg-n_elgqso)
# w = fracp * (n_elg-n_elgqso) / (n_elgqso - fracp * n_elgqso)
# TBD: get target files for main...
"""def get_coadd_weight(tileid):
    # parent sample
    #print(tileid,  glob("{}/survey/fiberassign/SV3/202?????/{:06d}-targ.fits".format(os.getenv("DESI_ROOT"), tileid)))
    fn = glob("{}/survey/fiberassign/SV3/202?????/{:06d}-targ.fits".format(os.getenv("DESI_ROOT"), tileid))[0]
    parent = fitsio.read(fn, columns = ["SV3_DESI_TARGET"])
    # assigned sample
    fn = glob("{}/survey/fiberassign/SV3/202?????/fiberassign-{:06d}.fits.gz".format(os.getenv("DESI_ROOT"), tileid))[0]
    assgn = fitsio.read(fn, columns = [dtkey])
    wdict = {}
    for tracer in elg_tracers:
        # parent
        elg = (parent[dtkey] & desi_mask[tracer]) > 0
        elgqso = (elg) & ((parent[dtkey] & desi_mask["QSO"]) > 0)
        fracp = elgqso.sum() / float(elg.sum())
        # assigned sample
        elg = (assgn[dtkey] & desi_mask[tracer]) > 0
        elgqso = (elg) & ((assgn[dtkey] & desi_mask["QSO"]) > 0)
        nelg, nelgqso = elg.sum(), elgqso.sum()
        #
        wdict[tracer] = np.round(fracp * (nelg - nelgqso) / (nelgqso - fracp*nelgqso),3)
    return wdict"""


def process_redux(tileid_thrunight_petal):
    #
    tileid, thrunight, petal = decode_tsp(tileid_thrunight_petal)
    print(tileid, thrunight, petal)
    # zbest filename
    zbfn = get_zbfn(tileid, thrunight, petal)
    # output filename for this zbest
    toutfn = get_outfn(tileid, thrunight, petal)
    # processing if zbest file present, and not already processed
    if os.path.isfile(zbfn) and not os.path.isfile(toutfn):
        #
        h = fits.open(zbfn)
        h["FIBERMAP"].data = h["FIBERMAP"].data[:500]
        ii_zb, ii_fm = match(h["ZBEST"].data["TARGETID"], h["FIBERMAP"].data["TARGETID"])
        if (len(ii_zb) != 500):
            sys.exit(zbfn+"; exiting")
        h["ZBEST"].data = h["ZBEST"].data[ii_zb]
        h["FIBERMAP"].data = h["FIBERMAP"].data[ii_fm]
        # tracer selection
        keep = (h["FIBERMAP"].data[dtkey] & desi_mask[args.tracer]) > 0
        h["ZBEST"].data = h["ZBEST"].data[keep]
        h["FIBERMAP"].data = h["FIBERMAP"].data[keep]
        ntracer = keep.sum()
        #
        mydict = {}
        # hascoadd?
        mydict["FN"] = [zbfn for i in range(ntracer)]
        mydict["THRUNIGHT"] = [thrunight for i in range(ntracer)]
        hascoadd = False
        if os.path.isfile(zbfn.replace("zbest","coadd")):
            tmph = fits.open(zbfn.replace("zbest","coadd"))
            extnames = [tmph[i].header["extname"] for i in range(1,len(tmph))]
            if "BRZ_WAVELENGTH" in extnames:
                hascoadd = True
            else:
                req_extnames = []
                for camera in ["B", "R", "Z"]:
                    req_extnames += [camera+"_WAVELENGTH", camera+"_FLUX", camera+"_IVAR"]
                if np.product([extname in extnames for extname in req_extnames]) == 1:
                    hascoadd = True
        mydict["HASCOADD"] = [hascoadd for i in range(ntracer)]
        # foii, foii_err
        if hascoadd:
            # mydict["FOII"], mydict["FOII_ERR"] = get_foii(zbfn.replace("zbest","coadd"), zbfn, targetids = h["ZBEST"].data["TARGETID"])
            scores = fits.open(zbfn.replace("zbest","coadd"))["SCORES"].data
            ii = np.array([np.where(scores["TARGETID"] == tid)[0][0] for tid in h["ZBEST"].data["TARGETID"]])
            if (scores["TARGETID"][ii] != h["ZBEST"].data["TARGETID"]).sum() > 0:
                sys.exit("issue with scores "+zbfn+"; exiting")
            for key in ["TSNR2_ELG", "TSNR2_LRG", "TSNR2_BGS"]:
                mydict[key] = scores[key][ii]
        else:
            # mydict["FOII"], mydict["FOII_ERR"] = -99 + np.zeros(ntracer), -99 + np.zeros(ntracer)
            for key in ["TSNR2_ELG", "TSNR2_LRG", "TSNR2_BGS"]:
                mydict[key] = -99 + np.zeros(ntracer)
        # # g, gr, rz, rw1, rw2
        # mydict["GTOT"], mydict["GFIB"], mydict["GR"], mydict["RZ"], mydict["RW1"], mydict["RW2"] = get_magcol(h["FIBERMAP"].data)
        # efftime
        efftime_spec, efftime_etc, efftime_gfa = get_efftime(zbfn)
        mydict["COADD_EFFTIME_SPEC"] = np.zeros(ntracer) + efftime_spec
        mydict["COADD_EFFTIME_ETC"] = np.zeros(ntracer) + efftime_etc
        mydict["COADD_EFFTIME_GFA"] = np.zeros(ntracer) + efftime_gfa
        # zdone (appears as a string in tiles.csv...)
        zdone = False
        keep = status_tiles["TILEID"] == tileid
        if keep.sum() > 0:
            zdone = status_tiles["ZDONE"][keep][0]
            zdone = zdone.lower() == "true"
        mydict["ZDONE"] = np.zeros(ntracer, dtype=bool) + zdone
        # weights for qso x elg/elg_hip/elg_lop
        # wdict = get_coadd_weight(tileid)
        # qso = (h["FIBERMAP"].data[dtkey] & desi_mask["QSO"]) > 0
        # for tracer in elg_tracers:
        #     mydict["COADD_WEIGHT_{}".format(tracer)] = np.ones(ntracer)
        #     mydict["COADD_WEIGHT_{}".format(tracer)][qso] = wdict[tracer]
        # writing to fits
        cols = []
        for extname in ["ZBEST", "FIBERMAP"]:
            cols += [h[extname].columns[key] for key in keys[extname]]
        for key, fmt in zip(okeys, ofmts):
            cols += [fits.Column(name=key, format=fmt, array=mydict[key])]
        h = fits.BinTableHDU.from_columns(fits.ColDefs(cols))

        print(toutfn)
        if not os.path.isdir(os.path.dirname(toutfn)):
            try:
                os.makedirs(os.path.dirname(toutfn))
            except:
                pass

        h.writeto(toutfn, overwrite=True)
    return True


#
tileid_thrunight_petals = []
for tileid, thrunight in zip(tileids, thrunights):
    tileid_thrunight_petals += [encode_tsp(tileid, thrunight, petal) for petal in range(10)]
print("len(tileid_thrunight_petals) = {}".format(len(tileid_thrunight_petals)))


# list of output files
tsp_outfns = []
for tileid_thrunight_petal in tileid_thrunight_petals:
    tileid, thrunight, petal = decode_tsp(tileid_thrunight_petal)
    tsp_outfns += [get_outfn(tileid, thrunight, petal)]

# cleaning first
#"""
for fn in tsp_outfns:
    if os.path.isfile(fn):
        os.remove(fn)
#"""

# AR parallel processing?
if args.numproc > 1:
    pool = sharedmem.MapReduce(np=args.numproc)
    with pool:
        _ = pool.map(process_redux, tileid_thrunight_petals)
if args.numproc == 1:
    for tileid_thrunight_petal in tileid_thrunight_petals:
        _ = process_redux(tileid_thrunight_petal)


# AR merging all processed tsp
hs = [Table.read(fn, format="fits") for fn in tsp_outfns if os.path.isfile(fn)]
h = vstack(hs)
# AR writing to fits
if os.path.isfile(outfits):
    os.remove(outfits)
fits.writeto(outfits, h.as_array(), overwrite=True)

