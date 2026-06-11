#!/usr/bin/env python
"""Build DR11 reference star catalogs for LRG/ELG pixel masks.

This is a DR11-configurable replacement for the DR9 scripts in:

  desi_mask/reference/
  desi_mask/lrg_mask/lrg_gaia_mask_v1.py
  desi_mask/elg_mask/gaia_elg_mask_v1.py
  desi_mask/lrg_mask/lrg_wise_mask_v1.py

The science choices are intentionally copied from the DR9 scripts:

  * Gaia EDR3 is trimmed to stars within 1.5 deg of a brick center.
  * Gaia mask_mag is min(G, predicted DECam z + 1), with G fallback.
  * Bright Tycho-2 stars missing from Gaia are added.
  * LRG/ELG Gaia radius-vs-magnitude relations are unchanged.
  * LRG WISE radius-vs-magnitude relation is unchanged.

Only the footprint/release-specific pieces are configurable.
"""

from __future__ import annotations

import argparse
import multiprocessing.pool
import os
from pathlib import Path
from typing import Iterable

import fitsio
import healpy as hp
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack, vstack


GAIA_EDR3_DIR = "/dvs_ro/cfs/cdirs/cosmo/data/gaia/edr3/healpix"
LEGACY_RELEASE_ROOT = "/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr11"

DEFAULT_USER = os.environ.get("USER", "rongpu")
DEFAULT_OUTPUT_DIR = (
    f"/global/cfs/cdirs/desi/users/{DEFAULT_USER}/desi_mask/dr11_reference_stars"
)

DEFAULT_TYCHO2_REFERENCE = (
    "/dvs_ro/cfs/cdirs/desi/users/rongpu/useful/tycho2-reference.fits"
)
DEFAULT_WISE_2MASS = (
    "/dvs_ro/cfs/cdirs/desi/users/rongpu/useful/"
    "unwise_bitmask_bright_stars/w1_bright-2mass.fits.gz"
)
DEFAULT_WISE_W1_13P3 = (
    "/dvs_ro/cfs/cdirs/desi/users/rongpu/useful/w1_bright-13.3.fits"
)
DEFAULT_GAIA_READ_PROCESSES = 128


STAGE_ORDER = [
    "trim_tycho2",
    "trim_gaia_g18",
    "trim_gaia_g14_pm",
    "predict_decam",
    "build_gaia_reference",
    "build_gaia_supplement",
    "build_lrg_gaia",
    "build_elg_gaia",
    "trim_wise_2mass",
    "trim_wise_faint",
    "combine_wise",
    "build_lrg_wise",
]


DECam_COEFFS = {
    "g": [
        -0.1125681175,
        0.3506376997,
        0.9082025788,
        -1.0078309266,
        -1.4212131445,
        4.5685722177,
        -4.5719415419,
        2.3816887292,
        -0.7162270722,
        0.1247021438,
        -0.0114938710,
        0.0003949585,
        0.0000051647,
    ],
    "r": [
        0.1431278873,
        -0.2999797766,
        -0.0553379742,
        0.1544273115,
        0.3068634689,
        -0.9499143903,
        0.9769739362,
        -0.4926704528,
        0.1272539574,
        -0.0133178183,
        -0.0008153813,
        0.0003094116,
        -0.0000198891,
    ],
    "z": [
        0.5173814296,
        -1.0450176704,
        0.1529797809,
        0.1856005222,
        -0.2366580132,
        0.1018331214,
        -0.0189673240,
        0.0012988354,
    ],
}


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Regenerate DR11 reference star catalogs for LRG/ELG masks."
    )
    parser.add_argument(
        "--release",
        default="dr11",
        help="Release label used in default output names and brick templates.",
    )
    parser.add_argument(
        "--release-root",
        default=LEGACY_RELEASE_ROOT,
        help="Legacy Survey release root.",
    )
    parser.add_argument(
        "--fields",
        nargs="+",
        default=["south", "north"],
        help="Photometric fields to include.",
    )
    parser.add_argument(
        "--brick-path-template",
        default="{release_root}/{field}/survey-bricks-{release}-{field}.fits.gz",
        help=(
            "Python format string for survey-brick files. Available keys are "
            "{release_root}, {release}, and {field}."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for all generated DR11 catalogs.",
    )
    parser.add_argument(
        "--gaia-dir",
        default=GAIA_EDR3_DIR,
        help="Gaia EDR3 healpix directory containing healpix-xxxxx.fits files.",
    )
    parser.add_argument(
        "--tycho2-reference",
        default=DEFAULT_TYCHO2_REFERENCE,
        help="All-sky Tycho-2 reference catalog produced by tycho2_reference.py.",
    )
    parser.add_argument(
        "--wise-2mass",
        default=DEFAULT_WISE_2MASS,
        help="AllWISE+2MASS bright-star catalog.",
    )
    parser.add_argument(
        "--wise-w1-13p3",
        default=DEFAULT_WISE_W1_13P3,
        help="AllWISE W1MPRO<13.3 catalog.",
    )
    parser.add_argument(
        "--old-gaia-mask",
        default=None,
        help=(
            "Optional old bright-star mask used to build a compatibility "
            "supplement. If omitted, the supplement stage is skipped."
        ),
    )
    parser.add_argument(
        "--stages",
        nargs="+",
        default=["all"],
        choices=["all"] + STAGE_ORDER,
        help="Stages to run.",
    )
    parser.add_argument(
        "--search-radius-deg",
        type=float,
        default=1.5,
        help="Distance from brick centers used to trim reference catalogs.",
    )
    parser.add_argument(
        "--healpix-preselect-radius-deg",
        type=float,
        default=4.0,
        help="Coarse healpix preselection radius for Gaia and faint WISE inputs.",
    )
    parser.add_argument(
        "--gaia-nside",
        type=int,
        default=32,
        help="Gaia source healpix nside.",
    )
    parser.add_argument(
        "--gaia-read-processes",
        type=positive_int,
        default=DEFAULT_GAIA_READ_PROCESSES,
        help=(
            "Number of worker processes for reading and trimming Gaia healpix "
            "files. Use 1 for serial reads."
        ),
    )
    parser.add_argument(
        "--wise-preselect-nside",
        type=int,
        default=32,
        help="Healpix nside for faint WISE footprint preselection.",
    )
    parser.add_argument(
        "--clobber",
        action="store_true",
        help="Overwrite existing outputs. Without this, existing stage outputs are kept.",
    )
    return parser.parse_args()


def output_paths(args: argparse.Namespace) -> dict[str, Path]:
    out = Path(args.output_dir)
    release = args.release
    return {
        "tycho2": out / f"tycho2-reference-{release}.fits",
        "gaia_g18": out / f"gaia_edr3_g_18_{release}.fits",
        "gaia_g14_pm": out / f"gaia_edr3_g_14_pm_{release}.fits",
        "gaia_decam": out / f"gaia_edr3_g_18_predict_decam_{release}.fits",
        "gaia_reference": out / f"gaia_reference_{release}.fits",
        "gaia_supplement": out / f"gaia_reference_suppl_{release}.fits",
        "lrg_gaia": out / f"gaia_lrg_mask_{release}_v1.fits",
        "elg_gaia": out / f"gaia_elg_mask_{release}_v1.fits",
        "wise_2mass": out / f"w1_bright-2mass-{release}.fits",
        "wise_faint": out / f"w1_bright-13.3-{release}.fits",
        "wise_combined": out / f"w1_bright-2mass-13.3-{release}.fits",
        "lrg_wise": out / f"w1_bright-2mass-lrg_mask_{release}_v1.fits",
    }


def maybe_skip(path: Path, clobber: bool) -> bool:
    if path.exists() and not clobber:
        print(f"{path} exists; skipping. Use --clobber to overwrite.")
        return True
    path.parent.mkdir(parents=True, exist_ok=True)
    return False


def write_table(table: Table, path: Path, clobber: bool) -> None:
    if maybe_skip(path, clobber):
        return
    print(f"Writing {len(table)} rows to {path}")
    table.write(path, overwrite=clobber)


def lower_column_names(table: Table) -> Table:
    table = table.copy()
    for col in list(table.colnames):
        lower = col.lower()
        if col != lower and lower not in table.colnames:
            table.rename_column(col, lower)
    return table


def get_col(table: Table, name: str) -> str:
    if name in table.colnames:
        return name
    upper = name.upper()
    if upper in table.colnames:
        return upper
    lower = name.lower()
    if lower in table.colnames:
        return lower
    raise KeyError(f"Column {name!r} not found in {table.colnames}")


def read_bricks(args: argparse.Namespace) -> Table:
    bricks = []
    columns = ["brickid", "brickname", "ra", "dec", "ra1", "ra2", "dec1", "dec2"]
    for field in args.fields:
        path = args.brick_path_template.format(
            release_root=args.release_root,
            release=args.release,
            field=field,
        )
        print(f"Reading bricks for {field}: {path}")
        brick = Table(fitsio.read(str(path)))
        brick = lower_column_names(brick)
        missing = [col for col in columns if col not in brick.colnames]
        if missing:
            raise ValueError(f"{path} is missing expected columns {missing}")
        brick["field"] = field
        bricks.append(brick[columns + ["field"]])

    bricks = vstack(bricks, join_type="exact")
    _, idx = np.unique(bricks["brickid"], return_index=True)
    bricks = bricks[np.sort(idx)]
    print(f"Unique DR11 bricks: {len(bricks)}")
    return bricks


def skycoord_from_table(table: Table, ra_col: str = "ra", dec_col: str = "dec") -> SkyCoord:
    return SkyCoord(table[ra_col] * u.degree, table[dec_col] * u.degree, frame="icrs")


def unique_sources_near_bricks(
    source_ra: np.ndarray,
    source_dec: np.ndarray,
    brick_sky: SkyCoord,
    radius_arcsec: float,
) -> np.ndarray:
    source_sky = SkyCoord(source_ra * u.degree, source_dec * u.degree, frame="icrs")
    idx_source, _, _, _ = brick_sky.search_around_sky(
        source_sky, seplimit=radius_arcsec * u.arcsec
    )
    return np.unique(idx_source)


def healpix_pixels_near_bricks(
    nside: int,
    brick_sky: SkyCoord,
    radius_arcsec: float,
) -> np.ndarray:
    npix = hp.nside2npix(nside)
    hp_ra, hp_dec = hp.pix2ang(nside, np.arange(npix), nest=True, lonlat=True)
    hp_sky = SkyCoord(hp_ra * u.degree, hp_dec * u.degree, frame="icrs")
    idx_hp, _, _, _ = brick_sky.search_around_sky(
        hp_sky, seplimit=radius_arcsec * u.arcsec
    )
    return np.unique(idx_hp)


def read_gaia_healpix(path: Path, rows: np.ndarray | None, columns: list[str]) -> Table:
    return Table(fitsio.read(str(path), rows=rows, columns=columns))


_GAIA_TRIM_WORKER_STATE: dict[str, object] = {}


def init_gaia_trim_worker(
    gaia_dir: str,
    brick_ra: np.ndarray,
    brick_dec: np.ndarray,
    search_radius: float,
    max_g_mag: float,
    first_pass_columns: list[str],
    output_columns: list[str],
) -> None:
    global _GAIA_TRIM_WORKER_STATE
    _GAIA_TRIM_WORKER_STATE = {
        "gaia_dir": gaia_dir,
        "brick_sky": SkyCoord(
            brick_ra * u.degree, brick_dec * u.degree, frame="icrs"
        ),
        "search_radius": search_radius,
        "max_g_mag": max_g_mag,
        "first_pass_columns": first_pass_columns,
        "output_columns": output_columns,
    }


def read_trimmed_gaia_healpix(
    task: tuple[int, int],
) -> tuple[int, int, Table | None]:
    state = _GAIA_TRIM_WORKER_STATE
    if not state:
        raise RuntimeError("Gaia trim worker is not initialized.")

    position, hp_index = task
    hp_index = int(hp_index)
    gaia_fn = Path(state["gaia_dir"]) / f"healpix-{hp_index:05d}.fits"

    tmp = read_gaia_healpix(
        gaia_fn, rows=None, columns=state["first_pass_columns"]
    )
    mag_mask = tmp["PHOT_G_MEAN_MAG"] < state["max_g_mag"]
    if np.sum(mag_mask) == 0:
        return position, hp_index, None

    idx = unique_sources_near_bricks(
        np.asarray(tmp["RA"][mag_mask]),
        np.asarray(tmp["DEC"][mag_mask]),
        state["brick_sky"],
        state["search_radius"],
    )
    if len(idx) == 0:
        return position, hp_index, None

    rows = np.where(mag_mask)[0][idx]
    table = read_gaia_healpix(
        gaia_fn, rows=rows, columns=state["output_columns"]
    )
    return position, hp_index, table


def read_trimmed_gaia_healpix_tables(
    gaia_pixels: np.ndarray,
    args: argparse.Namespace,
    bricks: Table,
    max_g_mag: float,
    first_pass_columns: list[str],
    output_columns: list[str],
    search_radius: float,
) -> list[Table]:
    tasks = [
        (position, int(hp_index)) for position, hp_index in enumerate(gaia_pixels)
    ]
    if len(tasks) == 0:
        return []

    worker_count = min(args.gaia_read_processes, len(tasks))
    init_args = (
        args.gaia_dir,
        np.asarray(bricks["ra"], dtype="f8"),
        np.asarray(bricks["dec"], dtype="f8"),
        search_radius,
        max_g_mag,
        first_pass_columns,
        output_columns,
    )
    tables_by_position: list[Table | None] = [None] * len(tasks)
    completed = 0

    print(f"Reading Gaia healpix files with {worker_count} process(es)")

    def store_result(result: tuple[int, int, Table | None]) -> None:
        nonlocal completed
        position, _, table = result
        completed += 1
        if completed == 1 or completed % 100 == 0 or completed == len(tasks):
            print(f"{completed}/{len(tasks)} Gaia healpix files processed")
        if table is not None:
            tables_by_position[position] = table

    if worker_count == 1:
        init_gaia_trim_worker(*init_args)
        for result in map(read_trimmed_gaia_healpix, tasks):
            store_result(result)
    else:
        with multiprocessing.pool.Pool(
            processes=worker_count,
            initializer=init_gaia_trim_worker,
            initargs=init_args,
        ) as pool:
            for result in pool.imap_unordered(
                read_trimmed_gaia_healpix, tasks, chunksize=1
            ):
                store_result(result)

    return [table for table in tables_by_position if table is not None]


def trim_gaia_catalog(
    args: argparse.Namespace,
    bricks: Table,
    max_g_mag: float,
    output_path: Path,
    include_pm: bool,
) -> None:
    if maybe_skip(output_path, args.clobber):
        return

    brick_sky = skycoord_from_table(bricks)
    search_radius = args.search_radius_deg * 3600.0
    init_radius = args.healpix_preselect_radius_deg * 3600.0

    print(f"Gaia healpix nside={args.gaia_nside}")
    print(f"Healpix pixel area: {hp.nside2pixarea(args.gaia_nside, degrees=True):.5f} deg2")
    gaia_pixels = healpix_pixels_near_bricks(args.gaia_nside, brick_sky, init_radius)
    print(f"Gaia healpix files near DR11 bricks: {len(gaia_pixels)}")

    first_pass_columns = ["RA", "DEC", "PHOT_G_MEAN_MAG"]
    output_columns = [
        "SOURCE_ID",
        "RA",
        "DEC",
        "PHOT_G_MEAN_MAG",
        "PHOT_BP_MEAN_MAG",
        "PHOT_RP_MEAN_MAG",
        "PHOT_G_MEAN_FLUX_OVER_ERROR",
        "ASTROMETRIC_EXCESS_NOISE",
    ]
    if include_pm:
        output_columns += [
            "PARALLAX",
            "PARALLAX_ERROR",
            "PMRA",
            "PMRA_ERROR",
            "PMDEC",
            "PMDEC_ERROR",
        ]

    gaia_tables = read_trimmed_gaia_healpix_tables(
        gaia_pixels=gaia_pixels,
        args=args,
        bricks=bricks,
        max_g_mag=max_g_mag,
        first_pass_columns=first_pass_columns,
        output_columns=output_columns,
        search_radius=search_radius,
    )

    if len(gaia_tables) == 0:
        raise RuntimeError("No Gaia sources found near DR11 bricks.")

    gaia = vstack(gaia_tables)
    print(f"Gaia rows before flux-over-error filter: {len(gaia)}")
    gaia = gaia[gaia["PHOT_G_MEAN_FLUX_OVER_ERROR"] > 0]
    print(f"Gaia rows after flux-over-error filter: {len(gaia)}")
    gaia.write(output_path, overwrite=args.clobber)


def trim_tycho2(args: argparse.Namespace, bricks: Table, paths: dict[str, Path]) -> None:
    output_path = paths["tycho2"]
    if maybe_skip(output_path, args.clobber):
        return

    tycho2 = Table(fitsio.read(args.tycho2_reference))
    tycho2["idx"] = np.arange(len(tycho2))

    brick_sky = skycoord_from_table(bricks)
    idx = unique_sources_near_bricks(
        np.asarray(tycho2["RA"]),
        np.asarray(tycho2["DEC"]),
        brick_sky,
        args.search_radius_deg * 3600.0,
    )
    print(f"Tycho-2 rows near DR11 bricks: {len(idx)} / {len(tycho2)}")
    tycho2 = tycho2[idx]
    tycho2.write(output_path, overwrite=args.clobber)


def trim_gaia_g18(args: argparse.Namespace, bricks: Table, paths: dict[str, Path]) -> None:
    trim_gaia_catalog(args, bricks, 18.0, paths["gaia_g18"], include_pm=False)


def trim_gaia_g14_pm(args: argparse.Namespace, bricks: Table, paths: dict[str, Path]) -> None:
    trim_gaia_catalog(args, bricks, 14.0, paths["gaia_g14_pm"], include_pm=True)


def predict_decam(args: argparse.Namespace, paths: dict[str, Path]) -> None:
    output_path = paths["gaia_decam"]
    if maybe_skip(output_path, args.clobber):
        return

    gaia = Table(fitsio.read(str(paths["gaia_g18"])))
    bp_rp = np.asarray(gaia["PHOT_BP_MEAN_MAG"] - gaia["PHOT_RP_MEAN_MAG"])
    bp_rp = np.clip(bp_rp, -0.5, 4.7)

    result = Table()
    invalid_color = (
        (gaia["PHOT_BP_MEAN_MAG"] == 0)
        | (gaia["PHOT_RP_MEAN_MAG"] == 0)
        | ~np.isfinite(bp_rp)
    )
    for band, coeffs in DECam_COEFFS.items():
        mag = np.array(gaia["PHOT_G_MEAN_MAG"], dtype="f8")
        for order, coeff in enumerate(coeffs):
            mag += coeff * bp_rp**order
        mag[invalid_color] = np.nan
        result[f"decam_mag_{band}"] = mag

    result.write(output_path, overwrite=args.clobber)


def compute_tycho_missing_from_gaia(tycho2: Table, gaia_pm: Table) -> Table:
    mask = tycho2["MAG_VT"] < 10.0
    mask &= (tycho2["MAG_VT"] != 0) | (tycho2["MAG_HP"] != 0)
    tycho2 = tycho2[mask].copy()
    print(f"Bright Tycho-2 candidates: {len(tycho2)}")

    valid_epoch = (tycho2["EPOCH_RA"] != 0) & (tycho2["EPOCH_DEC"] != 0)
    tycho_ra_j1991 = np.array(tycho2["RA"], dtype="f8")
    tycho_dec_j1991 = np.array(tycho2["DEC"], dtype="f8")
    cos_dec = np.cos(np.radians(np.asarray(tycho2["DEC"])))
    tycho_ra_j1991[valid_epoch] = (
        tycho2["RA"] - (tycho2["EPOCH_RA"] - 1991.0) * tycho2["PM_RA"] / 3600.0 / cos_dec
    )[valid_epoch]
    tycho_dec_j1991[valid_epoch] = (
        tycho2["DEC"] - (tycho2["EPOCH_DEC"] - 1991.0) * tycho2["PM_DEC"] / 3600.0
    )[valid_epoch]

    gaia_ra_j1991 = np.asarray(gaia_pm["RA"]) - (
        16.0
        * np.asarray(gaia_pm["PMRA"])
        * 1.0e-3
        / 3600.0
        / np.cos(np.radians(np.asarray(gaia_pm["DEC"])))
    )
    gaia_dec_j1991 = np.asarray(gaia_pm["DEC"]) - (
        16.0 * np.asarray(gaia_pm["PMDEC"]) * 1.0e-3 / 3600.0
    )

    tycho_sky = SkyCoord(
        tycho_ra_j1991 * u.degree, tycho_dec_j1991 * u.degree, frame="icrs"
    )
    gaia_sky = SkyCoord(gaia_ra_j1991 * u.degree, gaia_dec_j1991 * u.degree, frame="icrs")
    idx_tycho, _, _, _ = gaia_sky.search_around_sky(
        tycho_sky, seplimit=1.0 * u.arcsec
    )

    missing = np.ones(len(tycho2), dtype=bool)
    missing[np.unique(idx_tycho)] = False
    print(f"Tycho-2 stars missing from Gaia match: {np.sum(missing)} / {len(tycho2)}")
    tycho2 = tycho2[missing].copy()

    tycho2["mask_mag"] = np.minimum(tycho2["ggguess"], tycho2["zguess"] + 1.0)
    missing_z = np.isnan(tycho2["zguess"])
    use_vt = missing_z & (tycho2["MAG_VT"] != 0)
    use_hp = missing_z & (tycho2["MAG_VT"] == 0) & (tycho2["MAG_HP"] != 0)
    tycho2["mask_mag"][use_vt] = tycho2["MAG_VT"][use_vt]
    tycho2["mask_mag"][use_hp] = tycho2["MAG_HP"][use_hp]

    tycho2 = tycho2[["RA", "DEC", "mask_mag", "ggguess", "zguess"]]
    tycho2.rename_column("ggguess", "gaia_g")
    tycho2["is_tycho2"] = np.ones(len(tycho2), dtype=bool)
    return tycho2


def build_gaia_reference(args: argparse.Namespace, paths: dict[str, Path]) -> None:
    output_path = paths["gaia_reference"]
    if maybe_skip(output_path, args.clobber):
        return

    tycho2 = Table(fitsio.read(str(paths["tycho2"])))
    gaia_pm = Table(fitsio.read(str(paths["gaia_g14_pm"])))
    tycho_missing = compute_tycho_missing_from_gaia(tycho2, gaia_pm)

    gaia_columns = ["RA", "DEC", "PHOT_G_MEAN_MAG", "PHOT_G_MEAN_FLUX_OVER_ERROR"]
    gaia = Table(fitsio.read(str(paths["gaia_g18"]), columns=gaia_columns))
    decam = Table(fitsio.read(str(paths["gaia_decam"])))
    gaia = hstack([gaia, decam], join_type="exact")

    gaia["mask_mag"] = np.minimum(gaia["PHOT_G_MEAN_MAG"], gaia["decam_mag_z"] + 1.0)
    missing_z = np.isnan(gaia["decam_mag_z"])
    gaia["mask_mag"][missing_z] = gaia["PHOT_G_MEAN_MAG"][missing_z]
    print(f"Gaia rows with NaN mask_mag: {np.sum(np.isnan(gaia['mask_mag']))}")

    gaia = gaia[["RA", "DEC", "mask_mag", "PHOT_G_MEAN_MAG", "decam_mag_z"]]
    gaia.rename_columns(["PHOT_G_MEAN_MAG", "decam_mag_z"], ["gaia_g", "zguess"])
    gaia["is_tycho2"] = np.zeros(len(gaia), dtype=bool)

    ref = vstack([tycho_missing, gaia], join_type="exact")
    ref.write(output_path, overwrite=args.clobber)


def build_gaia_supplement(args: argparse.Namespace, paths: dict[str, Path]) -> None:
    if args.old_gaia_mask is None:
        print("No --old-gaia-mask supplied; skipping Gaia supplement.")
        return

    output_path = paths["gaia_supplement"]
    if maybe_skip(output_path, args.clobber):
        return

    old = Table(fitsio.read(args.old_gaia_mask, columns=["mask_mag"]))
    old_idx = np.where(old["mask_mag"] < 8.0)[0]
    old = Table(fitsio.read(args.old_gaia_mask, rows=old_idx))
    old = lower_column_names(old)

    new = Table(fitsio.read(str(paths["gaia_reference"]), columns=["RA", "DEC", "mask_mag"]))
    new = new[new["mask_mag"] < 10.0]

    old_ra = get_col(old, "ra")
    old_dec = get_col(old, "dec")
    old_sky = SkyCoord(old[old_ra] * u.degree, old[old_dec] * u.degree, frame="icrs")
    new_sky = SkyCoord(new["RA"] * u.degree, new["DEC"] * u.degree, frame="icrs")

    idx_new, idx_old, _, _ = old_sky.search_around_sky(
        new_sky, seplimit=5.0 * u.arcsec
    )
    old_unique = np.unique(idx_old)
    brightest_new = np.full(len(old_unique), np.inf)
    for ii, old_index in enumerate(old_unique):
        brightest_new[ii] = np.min(new["mask_mag"][idx_new[idx_old == old_index]])

    old_mag = np.asarray(old["mask_mag"][old_unique])
    keep = brightest_new - old_mag > 0.05
    supplement = old[old_unique[keep]]
    print(f"Gaia supplement rows: {len(supplement)}")
    supplement.write(output_path, overwrite=args.clobber)


def interp_log_radius(mags: Iterable[float], radii: Iterable[float]):
    mags = np.asarray(list(mags), dtype="f8")
    log_radii = np.log10(np.asarray(list(radii), dtype="f8"))

    def func(values: np.ndarray) -> np.ndarray:
        values = np.asarray(values, dtype="f8")
        clipped = np.interp(values, mags, log_radii)
        left = values < mags[0]
        right = values > mags[-1]
        if np.any(left):
            slope = (log_radii[1] - log_radii[0]) / (mags[1] - mags[0])
            clipped[left] = log_radii[0] + slope * (values[left] - mags[0])
        if np.any(right):
            slope = (log_radii[-1] - log_radii[-2]) / (mags[-1] - mags[-2])
            clipped[right] = log_radii[-1] + slope * (values[right] - mags[-1])
        return 10.0**clipped

    return func


def interp_linear_radius(mags: Iterable[float], radii: Iterable[float]):
    mags = np.asarray(list(mags), dtype="f8")
    radii = np.asarray(list(radii), dtype="f8")

    def func(values: np.ndarray) -> np.ndarray:
        values = np.asarray(values, dtype="f8")
        interp = np.interp(values, mags, radii)
        left = values < mags[0]
        right = values > mags[-1]
        if np.any(left):
            slope = (radii[1] - radii[0]) / (mags[1] - mags[0])
            interp[left] = radii[0] + slope * (values[left] - mags[0])
        if np.any(right):
            slope = (radii[-1] - radii[-2]) / (mags[-1] - mags[-2])
            interp[right] = radii[-1] + slope * (values[right] - mags[-1])
        return interp

    return func


def read_optional_gaia_supplement(path: Path) -> Table | None:
    if not path.exists():
        return None
    suppl = Table(fitsio.read(str(path)))
    suppl = lower_column_names(suppl)
    ra_col = get_col(suppl, "ra")
    dec_col = get_col(suppl, "dec")
    mag_col = get_col(suppl, "mask_mag")
    out = Table()
    out["RA"] = np.asarray(suppl[ra_col])
    out["DEC"] = np.asarray(suppl[dec_col])
    out["mask_mag"] = np.asarray(suppl[mag_col])
    return out


def read_gaia_reference_for_mask(paths: dict[str, Path]) -> Table:
    gaia = Table(fitsio.read(str(paths["gaia_reference"]), columns=["RA", "DEC", "mask_mag"]))
    suppl = read_optional_gaia_supplement(paths["gaia_supplement"])
    if suppl is not None:
        print(f"Adding Gaia supplement rows: {len(suppl)}")
        gaia = vstack([gaia, suppl], join_type="exact")
    return gaia


def build_lrg_gaia(args: argparse.Namespace, paths: dict[str, Path]) -> None:
    output_path = paths["lrg_gaia"]
    if maybe_skip(output_path, args.clobber):
        return

    gaia = read_gaia_reference_for_mask(paths)
    mask = gaia["mask_mag"] < 18.0
    gaia["radius_south"] = np.zeros(len(gaia), dtype="f8")
    gaia["radius_north"] = np.zeros(len(gaia), dtype="f8")

    f_south = interp_log_radius(
        [4.0, 9.0, 10.0, 10.5, 11.5, 12.0, 12.5, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 17.0, 18.0],
        [429.18637985, 80.95037032, 57.98737129, 36.80882682, 26.36735446, 25.29190318, 21.40616169, 15.33392671, 13.74150366, 13.56870306, 12.03092488, 11.10823009, 9.79334208, 7.01528803, 5.02527796],
    )
    f_north = interp_log_radius(
        [4.0, 9.0, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 17.0, 18.0],
        [429.18637985, 80.95037032, 60.0, 60.0, 60.0, 47.46123803, 38.68173428, 32.73883553, 27.70897871, 23.45188791, 19.84883862, 16.79934664, 13.67150555, 11.57107301, 7.83467367, 5.61223042, 4.02022236],
    )
    gaia["radius_south"][mask] = f_south(gaia["mask_mag"][mask])
    gaia["radius_north"][mask] = f_north(gaia["mask_mag"][mask])
    gaia.write(output_path, overwrite=args.clobber)


def build_elg_gaia(args: argparse.Namespace, paths: dict[str, Path]) -> None:
    output_path = paths["elg_gaia"]
    if maybe_skip(output_path, args.clobber):
        return

    gaia = read_gaia_reference_for_mask(paths)
    mask = gaia["mask_mag"] < 18.0
    gaia["radius_south"] = np.zeros(len(gaia), dtype="f8")
    gaia["radius_north"] = np.zeros(len(gaia), dtype="f8")

    mags = [4.0, 5.0, 6.0, 7.0] + np.arange(7.5, 18.5, 0.5).tolist()
    f_south = interp_log_radius(
        mags,
        [700.0, 500.0, 350.0, 260.0, 200.0, 180.0, 170.0, 120.0, 100.0, 75.0, 65.0, 55.0, 50.0, 40.0, 30.0, 25.0, 22.0, 19.0, 18.0, 14.0, 11.0, 10.0, 7.5, 6.0, 5.5, 4.2],
    )
    f_north = interp_log_radius(
        mags,
        [690.0, 590.0, 550.0, 510.0, 310.0, 270.0, 260.0, 250.0, 150.0, 120.0, 95.0, 85.0, 70.0, 65.0, 40.0, 35.0, 30.0, 25.0, 20.0, 18.0, 15.0, 12.0, 10.0, 10.0, 8.0, 7.0],
    )
    dr9_medium = lambda values: 1630.0 * 1.396 ** (-np.asarray(values))
    gaia["radius_south"][mask] = np.maximum(
        f_south(gaia["mask_mag"][mask]), dr9_medium(gaia["mask_mag"][mask])
    )
    gaia["radius_north"][mask] = np.maximum(
        f_north(gaia["mask_mag"][mask]), dr9_medium(gaia["mask_mag"][mask])
    )
    gaia.write(output_path, overwrite=args.clobber)


def trim_wise_2mass(args: argparse.Namespace, bricks: Table, paths: dict[str, Path]) -> None:
    output_path = paths["wise_2mass"]
    if maybe_skip(output_path, args.clobber):
        return

    wise = Table(fitsio.read(args.wise_2mass))
    brick_sky = skycoord_from_table(bricks)
    idx = unique_sources_near_bricks(
        np.asarray(wise["RA"]),
        np.asarray(wise["DEC"]),
        brick_sky,
        args.search_radius_deg * 3600.0,
    )
    print(f"WISE+2MASS rows near DR11 bricks: {len(idx)} / {len(wise)}")
    wise = wise[idx]
    wise.write(output_path, overwrite=args.clobber)


def trim_wise_faint(args: argparse.Namespace, bricks: Table, paths: dict[str, Path]) -> None:
    output_path = paths["wise_faint"]
    if maybe_skip(output_path, args.clobber):
        return

    wise = Table(fitsio.read(args.wise_w1_13p3))
    wise["idx"] = np.arange(len(wise))
    brick_sky = skycoord_from_table(bricks)

    init_radius = args.healpix_preselect_radius_deg * 3600.0
    hp_idx = healpix_pixels_near_bricks(args.wise_preselect_nside, brick_sky, init_radius)
    print(f"WISE preselection healpix pixels: {len(hp_idx)}")

    wise_hp = hp.ang2pix(
        args.wise_preselect_nside,
        np.asarray(wise["RA"]),
        np.asarray(wise["DEC"]),
        nest=True,
        lonlat=True,
    )
    wise = wise[np.isin(wise_hp, hp_idx)]
    print(f"WISE rows after healpix preselection: {len(wise)}")

    idx = unique_sources_near_bricks(
        np.asarray(wise["RA"]),
        np.asarray(wise["DEC"]),
        brick_sky,
        args.search_radius_deg * 3600.0,
    )
    wise = wise[idx]
    print(f"WISE rows near DR11 bricks: {len(wise)}")
    wise.write(output_path, overwrite=args.clobber)


def combine_wise(args: argparse.Namespace, paths: dict[str, Path]) -> None:
    output_path = paths["wise_combined"]
    if maybe_skip(output_path, args.clobber):
        return

    t1 = Table(fitsio.read(str(paths["wise_2mass"]), columns=["RA", "DEC", "W1MPRO"]))
    t2 = Table(fitsio.read(str(paths["wise_faint"]), columns=["RA", "DEC", "W1MPRO"]))
    t2 = t2[t2["W1MPRO"] >= 9.5]
    cat = vstack([t1, t2], join_type="exact")
    cat.write(output_path, overwrite=args.clobber)


def build_lrg_wise(args: argparse.Namespace, paths: dict[str, Path]) -> None:
    output_path = paths["lrg_wise"]
    if maybe_skip(output_path, args.clobber):
        return

    wise = Table(fitsio.read(str(paths["wise_combined"])))
    wise["w1ab"] = np.asarray(wise["W1MPRO"]) + 2.699
    wise["radius"] = np.zeros(len(wise), dtype="f8")

    f_radius = interp_linear_radius(
        [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0],
        [600.0, 600.0, 550.0, 500.0, 475.0, 425.0, 400.0, 400.0, 390.0, 392.5, 395.0, 370.0, 360.0, 330.0, 275.0, 240.0, 210.0, 165.0, 100.0, 75.0, 60.0],
    )
    mask = wise["w1ab"] < 10.0
    wise["radius"][mask] = f_radius(wise["w1ab"][mask])
    wise.write(output_path, overwrite=args.clobber)


def run_stage(
    stage: str,
    args: argparse.Namespace,
    bricks: Table,
    paths: dict[str, Path],
) -> None:
    print(f"\n=== {stage} ===")
    if stage == "trim_tycho2":
        trim_tycho2(args, bricks, paths)
    elif stage == "trim_gaia_g18":
        trim_gaia_g18(args, bricks, paths)
    elif stage == "trim_gaia_g14_pm":
        trim_gaia_g14_pm(args, bricks, paths)
    elif stage == "predict_decam":
        predict_decam(args, paths)
    elif stage == "build_gaia_reference":
        build_gaia_reference(args, paths)
    elif stage == "build_gaia_supplement":
        build_gaia_supplement(args, paths)
    elif stage == "build_lrg_gaia":
        build_lrg_gaia(args, paths)
    elif stage == "build_elg_gaia":
        build_elg_gaia(args, paths)
    elif stage == "trim_wise_2mass":
        trim_wise_2mass(args, bricks, paths)
    elif stage == "trim_wise_faint":
        trim_wise_faint(args, bricks, paths)
    elif stage == "combine_wise":
        combine_wise(args, paths)
    elif stage == "build_lrg_wise":
        build_lrg_wise(args, paths)
    else:
        raise ValueError(f"Unknown stage {stage}")


def main() -> None:
    args = parse_args()
    paths = output_paths(args)
    stages = STAGE_ORDER if "all" in args.stages else args.stages

    print("Output paths:")
    for key, path in paths.items():
        print(f"  {key}: {path}")

    bricks = read_bricks(args)
    for stage in stages:
        run_stage(stage, args, bricks, paths)


if __name__ == "__main__":
    main()

