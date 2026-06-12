#!/usr/bin/env python
"""Build DR11 per-brick LRG/ELG pixel masks.

This is the DR11 counterpart of the DR9 scripts in:

  desi_mask/lrg_mask/lrg_pixel_bitmask/gnu_parallel/
  desi_mask/elg_mask/elg_pixel_bitmask/
  desi_mask/all_tracers_custom_mask/

The DR9 scripts interpolated RA and Dec independently on each brick.  That is
fragile near RA wrap and the celestial poles, so this version interpolates unit
Cartesian vectors and renormalizes them before doing spherical distance tests.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import time
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path
from typing import Iterable

import fitsio
import numpy as np
from astropy import units as u
from astropy import wcs
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table, vstack
from scipy.interpolate import RectBivariateSpline
import multiprocessing


RELEASE_ROOT = "/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr11"
DEFAULT_USER = os.environ.get("USER", "rongpu")
DEFAULT_REFERENCE_DIR = (
    f"/global/cfs/cdirs/desi/users/{DEFAULT_USER}/desi_mask/dr11_reference_stars"
)
DEFAULT_OUTPUT_ROOT = (
    f"/global/cfs/cdirs/desi/users/{DEFAULT_USER}/desi_mask/dr11_pixel_masks"
)

STAGE_ORDER = [
    "all_custom",
    "elg_custom",
    "lrg_gaia",
    "lrg_wise",
    "elg_gaia",
    "combine_lrg",
    "combine_elg",
]

LRG_UNWISE_MASKBITS = [0, 1, 2, 3, 4, 6, 7]  # all except the HALO bit
TARGETING_MASKBITS = [1, 12, 13]

LRG_BITS = {
    "targeting": 0,
    "unwise": 1,
    "custom": 2,
    "gaia": 3,
    "wise": 4,
}

ELG_BITS = {
    "targeting": 0,
    "custom": 1,
    "gaia": 2,
    "elg_custom": 3,
}

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ALL_CUSTOM_MASK = "/global/cfs/cdirs/desicollab/users/rongpu/desi_mask/desi_custom_mask_v1.txt"
DEFAULT_ELG_CUSTOM_MASK = "/global/cfs/cdirs/desicollab/users/rongpu/desi_mask/elg_custom_mask_v1.txt"
DEBUG_BRICK_COUNT = 64

_WORKER_STATE: dict[str, object] = {}


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return parsed


def nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be >= 0")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build DR11 per-brick LRG/ELG pixel masks."
    )
    parser.add_argument("--release", default="dr11", help="Release label.")
    parser.add_argument(
        "--release-root",
        default=RELEASE_ROOT,
        help="Legacy Survey release root.",
    )
    parser.add_argument(
        "--fields",
        nargs="+",
        default=["south", "north"],
        help="Photometric fields to process.",
    )
    parser.add_argument(
        "--brick-path-template",
        default="{release_root}/{field}/survey-bricks-{release}-{field}.fits.gz",
        help=(
            "Format string for survey-brick files. Keys: release_root, release, field."
        ),
    )
    parser.add_argument(
        "--maskbits-path-template",
        default=(
            "{release_root}/{field}/coadd/{prefix}/{brickname}/"
            "legacysurvey-{brickname}-maskbits.fits.fz"
        ),
        help=(
            "Format string for per-brick maskbits files. Keys: release_root, "
            "release, field, prefix, brickname."
        ),
    )
    parser.add_argument(
        "--reference-dir",
        default=DEFAULT_REFERENCE_DIR,
        help="Directory containing DR11 reference-star catalogs.",
    )
    parser.add_argument(
        "--output-root",
        default=DEFAULT_OUTPUT_ROOT,
        help="Root directory for component and final pixel masks.",
    )
    parser.add_argument(
        "--stages",
        nargs="+",
        default=["all"],
        choices=["all"] + STAGE_ORDER,
        help="Stages to run.",
    )
    parser.add_argument(
        "--n-task",
        type=positive_int,
        default=1,
        help="Number of task splits for this stage.",
    )
    parser.add_argument(
        "--task-id",
        type=nonnegative_int,
        default=0,
        help="Task index in [0, n-task).",
    )
    parser.add_argument(
        "--processes",
        type=positive_int,
        default=128,
        help="Worker processes inside this task. Defaults to one Perlmutter CPU node.",
    )
    parser.add_argument(
        "--chunk-size",
        type=positive_int,
        default=400,
        help="Pixel chunk size for rasterization.",
    )
    parser.add_argument(
        "--spline-binsize",
        type=positive_int,
        default=1000,
        help="Coarse WCS sampling interval, in pixels, for Cartesian splines.",
    )
    parser.add_argument(
        "--star-batch-size",
        type=positive_int,
        default=256,
        help="Maximum stars per dot-product batch.",
    )
    parser.add_argument(
        "--shuffle-seed",
        type=int,
        default=4891,
        help="Seed used to shuffle bricks before task splitting.",
    )
    parser.add_argument(
        "--brickname",
        default=None,
        help="Optional single brick name for debugging.",
    )
    parser.add_argument(
        "--max-bricks",
        type=positive_int,
        default=None,
        help="Optional maximum number of bricks after filtering/splitting.",
    )
    parser.add_argument(
        "--debug-64-bricks",
        action="store_true",
        help="Debug mode: randomly select 64 total bricks before task splitting.",
    )
    parser.add_argument(
        "--clobber",
        action="store_true",
        help="Overwrite existing pixel-mask files.",
    )

    parser.add_argument(
        "--all-custom-mask-file",
        default=str(DEFAULT_ALL_CUSTOM_MASK),
        help="All-tracer custom mask text file.",
    )
    parser.add_argument(
        "--elg-custom-mask-file",
        default=str(DEFAULT_ELG_CUSTOM_MASK),
        help="ELG-specific custom mask text file.",
    )
    parser.add_argument(
        "--lrg-gaia-catalog",
        default=None,
        help="LRG Gaia mask reference catalog. Default is in --reference-dir.",
    )
    parser.add_argument(
        "--elg-gaia-catalog",
        default=None,
        help="ELG Gaia mask reference catalog. Default is in --reference-dir.",
    )
    parser.add_argument(
        "--lrg-wise-catalog",
        default=None,
        help="LRG WISE mask reference catalog. Default is in --reference-dir.",
    )

    parser.add_argument("--all-custom-dir", default=None)
    parser.add_argument("--elg-custom-dir", default=None)
    parser.add_argument("--lrg-gaia-dir", default=None)
    parser.add_argument("--lrg-wise-dir", default=None)
    parser.add_argument("--elg-gaia-dir", default=None)
    parser.add_argument("--lrg-output-dir", default=None)
    parser.add_argument("--elg-output-dir", default=None)

    args = parser.parse_args()
    if args.task_id >= args.n_task:
        parser.error("--task-id must be smaller than --n-task")
    return fill_default_paths(args)


def fill_default_paths(args: argparse.Namespace) -> argparse.Namespace:
    reference_dir = Path(args.reference_dir)
    output_root = Path(args.output_root)
    release = args.release

    if args.lrg_gaia_catalog is None:
        args.lrg_gaia_catalog = str(reference_dir / f"gaia_lrg_mask_{release}_v1.fits")
    if args.elg_gaia_catalog is None:
        args.elg_gaia_catalog = str(reference_dir / f"gaia_elg_mask_{release}_v1.fits")
    if args.lrg_wise_catalog is None:
        args.lrg_wise_catalog = str(
            reference_dir / f"w1_bright-2mass-lrg_mask_{release}_v1.fits"
        )

    args.all_custom_dir = args.all_custom_dir or str(
        output_root / "all_tracers_custom_mask" / "v1"
    )
    args.elg_custom_dir = args.elg_custom_dir or str(
        output_root / "elg" / "dev" / "elg_custom" / "v1"
    )
    args.lrg_gaia_dir = args.lrg_gaia_dir or str(
        output_root / "lrg" / "dev" / "gaiamask" / "v1"
    )
    args.lrg_wise_dir = args.lrg_wise_dir or str(
        output_root / "lrg" / "dev" / "wisemask" / "v1"
    )
    args.elg_gaia_dir = args.elg_gaia_dir or str(
        output_root / "elg" / "dev" / "gaiamask" / "v1"
    )
    args.lrg_output_dir = args.lrg_output_dir or str(output_root / "lrg" / "v1")
    args.elg_output_dir = args.elg_output_dir or str(output_root / "elg" / "v1")
    return args


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


def read_fits_columns(path: str, columns: Iterable[str]) -> Table:
    with fitsio.FITS(path) as hdul:
        available = list(hdul[1].get_colnames())
        lookup = {col.lower(): col for col in available}
        selected = []
        for col in columns:
            key = col.lower()
            if key not in lookup:
                raise KeyError(f"Column {col!r} not found in {available}")
            selected.append(lookup[key])
        return Table(hdul[1].read(columns=selected))


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
        brick = Table(fitsio.read(path))
        brick = lower_column_names(brick)
        missing = [col for col in columns if col not in brick.colnames]
        if missing:
            raise ValueError(f"{path} is missing expected columns {missing}")
        brick["field"] = field
        bricks.append(brick[columns + ["field"]])

    bricks = vstack(bricks, join_type="exact")
    if args.brickname is not None:
        bricks = bricks[np.asarray(bricks["brickname"]).astype(str) == args.brickname]
        if len(bricks) == 0:
            raise ValueError(f"Brick {args.brickname!r} was not found.")
    return bricks


def unique_bricks(bricks: Table) -> Table:
    _, idx = np.unique(bricks["brickid"], return_index=True)
    return bricks[np.sort(idx)]


def select_debug_bricks(bricks: Table, args: argparse.Namespace) -> Table:
    if not args.debug_64_bricks:
        return bricks

    sample_size = min(DEBUG_BRICK_COUNT, len(bricks))
    if sample_size == len(bricks):
        selected = np.arange(len(bricks))
    else:
        rng = np.random.default_rng(args.shuffle_seed)
        selected = rng.permutation(len(bricks))[:sample_size]
    print(f"Debug brick sample: {sample_size}/{len(bricks)} bricks")
    return bricks[selected]


def select_task_bricks(bricks: Table, args: argparse.Namespace) -> Table:
    rng = np.random.default_rng(args.shuffle_seed)
    order = rng.permutation(len(bricks))
    chunks = np.array_split(order, args.n_task)
    selected = chunks[args.task_id]
    if args.max_bricks is not None:
        selected = selected[: args.max_bricks]
    return bricks[selected]


def maskbits_path(args: argparse.Namespace, brick) -> Path:
    brickname = str(brick["brickname"])
    return Path(
        args.maskbits_path_template.format(
            release_root=args.release_root,
            release=args.release,
            field=str(brick["field"]),
            prefix=brickname[:3],
            brickname=brickname,
        )
    )


def all_custom_path(args: argparse.Namespace, brick) -> Path:
    brickname = str(brick["brickname"])
    return (
        Path(args.all_custom_dir)
        / "coadd"
        / brickname[:3]
        / brickname
        / f"{brickname}-custommask.fits.gz"
    )


def elg_custom_path(args: argparse.Namespace, brick) -> Path:
    brickname = str(brick["brickname"])
    return (
        Path(args.elg_custom_dir)
        / "coadd"
        / brickname[:3]
        / brickname
        / f"{brickname}-elg-custommask.fits.gz"
    )


def field_component_path(root: str, suffix: str, brick) -> Path:
    brickname = str(brick["brickname"])
    return (
        Path(root)
        / str(brick["field"])
        / "coadd"
        / brickname[:3]
        / brickname
        / f"{brickname}-{suffix}.fits.gz"
    )


def lrg_gaia_path(args: argparse.Namespace, brick) -> Path:
    return field_component_path(args.lrg_gaia_dir, "gaiamask", brick)


def lrg_wise_path(args: argparse.Namespace, brick) -> Path:
    return field_component_path(args.lrg_wise_dir, "wisemask", brick)


def elg_gaia_path(args: argparse.Namespace, brick) -> Path:
    return field_component_path(args.elg_gaia_dir, "gaiamask", brick)


def lrg_output_path(args: argparse.Namespace, brick) -> Path:
    return field_component_path(args.lrg_output_dir, "lrgmask", brick)


def elg_output_path(args: argparse.Namespace, brick) -> Path:
    return field_component_path(args.elg_output_dir, "elgmask", brick)


def radec_to_unit(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    cos_dec = np.cos(dec_rad)
    return np.stack(
        [
            np.cos(ra_rad) * cos_dec,
            np.sin(ra_rad) * cos_dec,
            np.sin(dec_rad),
        ],
        axis=0,
    )


def normalize_vectors(vec: np.ndarray) -> np.ndarray:
    norm = np.sqrt(np.sum(vec * vec, axis=0))
    return vec / norm


def unit_to_radec(vec: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    vec = normalize_vectors(vec)
    ra = (np.degrees(np.arctan2(vec[1], vec[0])) + 360.0) % 360.0
    dec = np.degrees(np.arcsin(np.clip(vec[2], -1.0, 1.0)))
    return ra, dec


def angular_separation_arcsec(vec1: np.ndarray, vec2: np.ndarray) -> np.ndarray:
    dot = np.sum(vec1 * vec2, axis=0)
    dot = np.clip(dot, -1.0, 1.0)
    return np.degrees(np.arccos(dot)) * 3600.0


def read_wcs_header(path: Path) -> tuple[dict[str, object], wcs.WCS, int, int]:
    with fits.open(path, memmap=False) as hdul:
        header = hdul[1].header.copy()

    header_keywords = [
        "CTYPE1",
        "CTYPE2",
        "CRVAL1",
        "CRVAL2",
        "CRPIX1",
        "CRPIX2",
        "CD1_1",
        "CD1_2",
        "CD2_1",
        "CD2_2",
    ]
    hdict = {keyword: header[keyword] for keyword in header_keywords}
    if hdict["CTYPE1"] != "RA---TAN" or hdict["CTYPE2"] != "DEC--TAN":
        raise ValueError(f"{path} does not use RA---TAN/DEC--TAN WCS.")

    naxis1 = int(header["NAXIS1"])
    naxis2 = int(header["NAXIS2"])
    return hdict, wcs.WCS(header), naxis1, naxis2


def build_cartesian_interpolators(
    ww: wcs.WCS, naxis1: int, naxis2: int, binsize: int
) -> list[RectBivariateSpline]:
    x_spline = np.arange(-binsize, naxis1 + 2 * binsize, binsize, dtype="f8")
    y_spline = np.arange(-binsize, naxis2 + 2 * binsize, binsize, dtype="f8")
    xx, yy = np.meshgrid(x_spline, y_spline)
    ra_grid, dec_grid = ww.wcs_pix2world(xx, yy, 0)
    vec_grid = radec_to_unit(ra_grid, dec_grid)
    return [
        RectBivariateSpline(y_spline, x_spline, vec_grid[ii])
        for ii in range(3)
    ]


def eval_cartesian_interpolators(
    interpolators: list[RectBivariateSpline], y_pix: np.ndarray, x_pix: np.ndarray
) -> np.ndarray:
    vec = np.stack([interp(y_pix, x_pix) for interp in interpolators], axis=0)
    return normalize_vectors(vec)


def wcs_image_center_and_radius(
    ww: wcs.WCS, naxis1: int, naxis2: int
) -> tuple[np.ndarray, float]:
    center_x = (naxis1 - 1) / 2.0
    center_y = (naxis2 - 1) / 2.0
    center_ra, center_dec = ww.wcs_pix2world([center_x], [center_y], 0)
    center_vec = radec_to_unit(np.asarray(center_ra), np.asarray(center_dec))

    x = np.array([0.0, 0.0, naxis1 - 1.0, naxis1 - 1.0, center_x, center_x, 0.0, naxis1 - 1.0])
    y = np.array([0.0, naxis2 - 1.0, 0.0, naxis2 - 1.0, 0.0, naxis2 - 1.0, center_y, center_y])
    ra, dec = ww.wcs_pix2world(x, y, 0)
    edge_vec = radec_to_unit(np.asarray(ra), np.asarray(dec))
    radius = float(np.max(angular_separation_arcsec(center_vec, edge_vec)))
    return center_vec[:, 0], radius + 1.0


def single_skycoord(ra: float, dec: float) -> SkyCoord:
    return SkyCoord([ra] * u.degree, [dec] * u.degree, frame="icrs")


def skycoord_from_arrays(ra: np.ndarray, dec: np.ndarray) -> SkyCoord | None:
    if len(ra) == 0:
        return None
    return SkyCoord(ra * u.degree, dec * u.degree, frame="icrs")


def brick_table_center_and_radius(brick) -> tuple[SkyCoord, float]:
    center_scalar = SkyCoord(float(brick["ra"]) * u.degree, float(brick["dec"]) * u.degree)
    corners_ra = np.array(
        [brick["ra1"], brick["ra1"], brick["ra2"], brick["ra2"]], dtype="f8"
    )
    corners_dec = np.array(
        [brick["dec1"], brick["dec2"], brick["dec1"], brick["dec2"]], dtype="f8"
    )
    corners = SkyCoord(corners_ra * u.degree, corners_dec * u.degree)
    radius = float(np.max(center_scalar.separation(corners).to_value(u.arcsec)))
    return single_skycoord(float(brick["ra"]), float(brick["dec"])), radius + 1.0


def read_star_catalog(path: str, radius_column: str | None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    columns = ["RA", "DEC"]
    if radius_column is not None:
        columns.append(radius_column)
    table = read_fits_columns(path, columns)
    ra_col = get_col(table, "RA")
    dec_col = get_col(table, "DEC")
    if radius_column is None:
        radius = np.zeros(len(table), dtype="f8")
    else:
        radius_col = get_col(table, radius_column)
        radius = np.asarray(table[radius_col], dtype="f8")
    ra = np.asarray(table[ra_col], dtype="f8")
    dec = np.asarray(table[dec_col], dtype="f8")
    keep = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(radius) & (radius > 0)
    return ra[keep], dec[keep], radius[keep]


def parse_custom_mask(path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    circular: list[list[float]] = []
    rectangular: list[list[float]] = []
    with open(path, "r", encoding="utf-8") as mask_file:
        for raw_line in mask_file:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            values = [float(value) for value in line.split(",")]
            if len(values) == 3:
                circular.append(values)
            elif len(values) == 4:
                rectangular.append(values)
            else:
                raise ValueError(f"Could not parse custom-mask line: {raw_line!r}")

    if circular:
        circ = np.asarray(circular, dtype="f8")
        ra, dec, radius = circ[:, 0], circ[:, 1], circ[:, 2]
        keep = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(radius) & (radius > 0)
        ra, dec, radius = ra[keep], dec[keep], radius[keep]
    else:
        ra = dec = radius = np.asarray([], dtype="f8")

    if rectangular:
        rectangles = np.asarray(rectangular, dtype="f8")
    else:
        rectangles = np.empty((0, 4), dtype="f8")
    return ra, dec, radius, rectangles


def rectangle_overlaps_brick(rectangle: np.ndarray, brick) -> bool:
    ramin, ramax, decmin, decmax = rectangle
    if decmax <= float(brick["dec1"]) or decmin >= float(brick["dec2"]):
        return False
    brick_ra1 = float(brick["ra1"])
    brick_ra2 = float(brick["ra2"])
    if ramin <= ramax:
        return ramin < brick_ra2 and ramax > brick_ra1
    return (ramin < brick_ra2 and 360.0 > brick_ra1) or (
        0.0 < brick_ra2 and ramax > brick_ra1
    )


def rectangle_mask(ra: np.ndarray, dec: np.ndarray, rectangles: np.ndarray) -> np.ndarray:
    mask = np.zeros(len(ra), dtype=bool)
    for ramin, ramax, decmin, decmax in rectangles:
        in_dec = (dec > decmin) & (dec < decmax)
        if ramin <= ramax:
            in_ra = (ra > ramin) & (ra < ramax)
        else:
            in_ra = (ra > ramin) | (ra < ramax)
        mask |= in_ra & in_dec
    return mask


def init_worker(state: dict[str, object]) -> None:
    global _WORKER_STATE
    _WORKER_STATE = state


def has_circular_overlap(brick, sky: SkyCoord | None, radius: np.ndarray) -> bool:
    if sky is None or len(radius) == 0:
        return False
    center, brick_radius = brick_table_center_and_radius(brick)
    max_radius = float(np.max(radius))
    _, idx, d2d, _ = sky.search_around_sky(
        center, seplimit=(max_radius + brick_radius) * u.arcsec
    )
    if len(idx) == 0:
        return False
    sep = np.asarray(d2d.to_value(u.arcsec))
    return bool(np.any(sep < radius[idx] + brick_radius))


def rasterize_circular_and_rectangular_masks(
    brick,
    header_path: Path,
    output_path: Path,
    args: argparse.Namespace,
    star_ra: np.ndarray,
    star_dec: np.ndarray,
    star_radius: np.ndarray,
    rectangles: np.ndarray,
    star_sky: SkyCoord | None,
    circle_bit: int | None,
    rectangle_bit: int | None,
    skip_empty: bool,
    precheck_overlap: bool,
) -> str:
    if output_path.exists() and not args.clobber:
        return "exists"

    if star_sky is None and len(star_ra) > 0:
        star_sky = skycoord_from_arrays(star_ra, star_dec)

    if precheck_overlap:
        rect_overlap = any(rectangle_overlaps_brick(rect, brick) for rect in rectangles)
        circ_overlap = has_circular_overlap(brick, star_sky, star_radius)
        if not rect_overlap and not circ_overlap:
            return "empty-skip"

    hdict, ww, naxis1, naxis2 = read_wcs_header(header_path)
    interpolators = build_cartesian_interpolators(
        ww, naxis1, naxis2, args.spline_binsize
    )
    center_vec, brick_radius = wcs_image_center_and_radius(ww, naxis1, naxis2)

    bitmask = np.zeros((naxis2, naxis1), dtype=np.uint8)
    brick_indices = np.asarray([], dtype=int)
    if star_sky is not None:
        center_ra, center_dec = unit_to_radec(center_vec[:, None])
        center_sky = single_skycoord(float(center_ra[0]), float(center_dec[0]))
        max_radius = float(np.max(star_radius))
        _, idx, d2d, _ = star_sky.search_around_sky(
            center_sky, seplimit=(max_radius + brick_radius) * u.arcsec
        )
        if len(idx) > 0:
            sep = np.asarray(d2d.to_value(u.arcsec))
            keep = sep < star_radius[idx] + brick_radius
            brick_indices = idx[keep]

    brick_star_ra = star_ra[brick_indices]
    brick_star_dec = star_dec[brick_indices]
    brick_star_radius = star_radius[brick_indices]
    brick_star_vec = radec_to_unit(brick_star_ra, brick_star_dec)
    brick_star_sky = None
    if len(brick_star_ra) > 0:
        brick_star_sky = SkyCoord(
            brick_star_ra * u.degree, brick_star_dec * u.degree, frame="icrs"
        )
        max_brick_star_radius = float(np.max(brick_star_radius))
    else:
        max_brick_star_radius = 0.0

    for y0 in range(0, naxis2, args.chunk_size):
        y1 = min(y0 + args.chunk_size, naxis2)
        y_pix = np.arange(y0, y1, dtype="f8")
        for x0 in range(0, naxis1, args.chunk_size):
            x1 = min(x0 + args.chunk_size, naxis1)
            x_pix = np.arange(x0, x1, dtype="f8")
            pix_vec = eval_cartesian_interpolators(interpolators, y_pix, x_pix)
            flat_vec = pix_vec.reshape(3, -1)
            chunk = bitmask[y0:y1, x0:x1]

            if rectangle_bit is not None and len(rectangles) > 0:
                ra, dec = unit_to_radec(flat_vec)
                rect = rectangle_mask(ra, dec, rectangles).reshape(y1 - y0, x1 - x0)
                chunk[rect] |= np.uint8(1 << rectangle_bit)

            if circle_bit is None or brick_star_sky is None:
                continue

            center = normalize_vectors(np.mean(flat_vec, axis=1)[:, None])[:, 0]
            chunk_radius = float(
                np.max(angular_separation_arcsec(center[:, None], flat_vec))
            ) + 1.0
            center_ra, center_dec = unit_to_radec(center[:, None])
            center_sky = single_skycoord(float(center_ra[0]), float(center_dec[0]))

            _, idx, d2d, _ = brick_star_sky.search_around_sky(
                center_sky,
                seplimit=(max_brick_star_radius + chunk_radius) * u.arcsec,
            )
            if len(idx) == 0:
                continue
            sep = np.asarray(d2d.to_value(u.arcsec))
            keep = sep < brick_star_radius[idx] + chunk_radius
            idx = idx[keep]
            if len(idx) == 0:
                continue

            circle = np.zeros(flat_vec.shape[1], dtype=bool)
            for start in range(0, len(idx), args.star_batch_size):
                sub_idx = idx[start : start + args.star_batch_size]
                star_vec = brick_star_vec[:, sub_idx]
                cos_radius = np.cos(np.radians(brick_star_radius[sub_idx] / 3600.0))
                dot = flat_vec.T @ star_vec
                circle |= np.any(dot > cos_radius[None, :], axis=1)
            circle = circle.reshape(y1 - y0, x1 - x0)
            chunk[circle] |= np.uint8(1 << circle_bit)

    if skip_empty and np.all(bitmask == 0):
        return "empty-skip"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fitsio.write(str(output_path), bitmask, compress="GZIP", header=hdict, clobber=True)
    return "wrote"


def worker_component(brick_index: int) -> tuple[str, str]:
    state = _WORKER_STATE
    args = state["args"]
    bricks = state["bricks"]
    brick = bricks[brick_index]
    stage = state["stage"]

    if stage == "all_custom":
        output_path = all_custom_path(args, brick)
        header_path = maskbits_path(args, brick)
        status = rasterize_circular_and_rectangular_masks(
            brick=brick,
            header_path=header_path,
            output_path=output_path,
            args=args,
            star_ra=state["star_ra"],
            star_dec=state["star_dec"],
            star_radius=state["star_radius"],
            rectangles=state["rectangles"],
            star_sky=state["star_sky"],
            circle_bit=0,
            rectangle_bit=1,
            skip_empty=True,
            precheck_overlap=True,
        )
    elif stage == "elg_custom":
        output_path = elg_custom_path(args, brick)
        header_path = maskbits_path(args, brick)
        status = rasterize_circular_and_rectangular_masks(
            brick=brick,
            header_path=header_path,
            output_path=output_path,
            args=args,
            star_ra=state["star_ra"],
            star_dec=state["star_dec"],
            star_radius=state["star_radius"],
            rectangles=state["rectangles"],
            star_sky=state["star_sky"],
            circle_bit=0,
            rectangle_bit=0,
            skip_empty=True,
            precheck_overlap=True,
        )
    elif stage == "lrg_gaia":
        output_path = lrg_gaia_path(args, brick)
        header_path = maskbits_path(args, brick)
        status = rasterize_circular_and_rectangular_masks(
            brick=brick,
            header_path=header_path,
            output_path=output_path,
            args=args,
            star_ra=state["star_ra"],
            star_dec=state["star_dec"],
            star_radius=state["star_radius"],
            rectangles=np.empty((0, 4), dtype="f8"),
            star_sky=state["star_sky"],
            circle_bit=0,
            rectangle_bit=None,
            skip_empty=False,
            precheck_overlap=False,
        )
    elif stage == "lrg_wise":
        output_path = lrg_wise_path(args, brick)
        header_path = maskbits_path(args, brick)
        status = rasterize_circular_and_rectangular_masks(
            brick=brick,
            header_path=header_path,
            output_path=output_path,
            args=args,
            star_ra=state["star_ra"],
            star_dec=state["star_dec"],
            star_radius=state["star_radius"],
            rectangles=np.empty((0, 4), dtype="f8"),
            star_sky=state["star_sky"],
            circle_bit=0,
            rectangle_bit=None,
            skip_empty=False,
            precheck_overlap=False,
        )
    elif stage == "elg_gaia":
        output_path = elg_gaia_path(args, brick)
        header_path = maskbits_path(args, brick)
        status = rasterize_circular_and_rectangular_masks(
            brick=brick,
            header_path=header_path,
            output_path=output_path,
            args=args,
            star_ra=state["star_ra"],
            star_dec=state["star_dec"],
            star_radius=state["star_radius"],
            rectangles=np.empty((0, 4), dtype="f8"),
            star_sky=state["star_sky"],
            circle_bit=0,
            rectangle_bit=None,
            skip_empty=False,
            precheck_overlap=False,
        )
    else:
        raise ValueError(f"Unsupported component stage {stage}")

    return str(brick["brickname"]), status


def add_maskbits(bitmask: np.ndarray, source: np.ndarray, bits: Iterable[int], out_bit: int) -> None:
    mask = np.zeros(source.shape, dtype=bool)
    for bit in bits:
        mask |= (source & (1 << bit)) > 0
    bitmask[mask] |= np.uint8(1 << out_bit)


def read_optional_component(path: Path) -> np.ndarray | None:
    if not path.exists():
        return None
    return fitsio.read(str(path)).astype(np.uint8)


def worker_combine(brick_index: int) -> tuple[str, str]:
    state = _WORKER_STATE
    args = state["args"]
    bricks = state["bricks"]
    stage = state["stage"]
    brick = bricks[brick_index]

    brickname = str(brick["brickname"])
    if stage == "combine_lrg":
        output_path = lrg_output_path(args, brick)
    elif stage == "combine_elg":
        output_path = elg_output_path(args, brick)
    else:
        raise ValueError(f"Unsupported combine stage {stage}")

    if output_path.exists() and not args.clobber:
        return brickname, "exists"

    input_path = maskbits_path(args, brick)
    hdict, _, _, _ = read_wcs_header(input_path)
    maskbits = fitsio.read(str(input_path), ext="MASKBITS")
    bitmask = np.zeros(maskbits.shape, dtype=np.uint8)

    add_maskbits(bitmask, maskbits, TARGETING_MASKBITS, 0)

    if stage == "combine_lrg":
        wisem1 = fitsio.read(str(input_path), ext="WISEM1")
        add_maskbits(bitmask, wisem1, LRG_UNWISE_MASKBITS, LRG_BITS["unwise"])

        custom = read_optional_component(all_custom_path(args, brick))
        if custom is not None:
            bitmask[custom != 0] |= np.uint8(1 << LRG_BITS["custom"])

        gaia = fitsio.read(str(lrg_gaia_path(args, brick))).astype(np.uint8)
        bitmask[gaia != 0] |= np.uint8(1 << LRG_BITS["gaia"])

        wise = fitsio.read(str(lrg_wise_path(args, brick))).astype(np.uint8)
        bitmask[wise != 0] |= np.uint8(1 << LRG_BITS["wise"])
    else:
        custom = read_optional_component(all_custom_path(args, brick))
        if custom is not None:
            bitmask[custom != 0] |= np.uint8(1 << ELG_BITS["custom"])

        elg_custom = read_optional_component(elg_custom_path(args, brick))
        if elg_custom is not None:
            bitmask[elg_custom != 0] |= np.uint8(1 << ELG_BITS["elg_custom"])

        gaia = fitsio.read(str(elg_gaia_path(args, brick))).astype(np.uint8)
        bitmask[gaia != 0] |= np.uint8(1 << ELG_BITS["gaia"])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fitsio.write(str(output_path), bitmask, compress="GZIP", header=hdict, clobber=True)
    return brickname, "wrote"


def shutdown_executor(executor, wait: bool) -> None:
    try:
        executor.shutdown(wait=wait, cancel_futures=not wait)
    except TypeError:
        executor.shutdown(wait=wait)


def process_bricks(
    stage: str,
    args: argparse.Namespace,
    bricks: Table,
    worker,
    state: dict[str, object],
) -> None:
    selected = select_task_bricks(bricks, args)
    processes = min(args.processes, len(selected))
    print(
        f"{stage}: task {args.task_id}/{args.n_task}; "
        f"{len(selected)} bricks; processes={processes}"
    )
    state = dict(state)
    state.update({"stage": stage, "args": args, "bricks": selected})

    time_start = time.time()
    counts: dict[str, int] = {}

    def store(result: tuple[str, str]) -> None:
        _, status = result
        counts[status] = counts.get(status, 0) + 1
        done = sum(counts.values())
        if done == 1 or done % 1000 == 0 or done == len(selected):
            print(f"{stage}: {done}/{len(selected)} {counts}", flush=True)

    if len(selected) == 0:
        print(f"{stage}: no bricks selected")
        return

    if processes == 1:
        init_worker(state)
        for index in range(len(selected)):
            store(worker(index))
    else:
        if multiprocessing is None:
            raise RuntimeError("multiprocessing is unavailable")
        executor = concurrent.futures.ProcessPoolExecutor(
            max_workers=processes,
            initializer=init_worker,
            initargs=(state,),
        )
        try:
            next_index = 0
            futures = {}

            def submit_next() -> None:
                nonlocal next_index
                if next_index < len(selected):
                    futures[executor.submit(worker, next_index)] = next_index
                    next_index += 1

            max_pending = min(len(selected), processes * 2)
            for _ in range(max_pending):
                submit_next()

            while futures:
                done, _ = concurrent.futures.wait(
                    futures, return_when=concurrent.futures.FIRST_COMPLETED
                )
                for future in done:
                    index = futures.pop(future)
                    try:
                        store(future.result())
                    except BrokenProcessPool:
                        raise
                    except Exception as exc:
                        brickname = str(selected[index]["brickname"])
                        raise RuntimeError(
                            f"{stage}: worker failed while processing brick {brickname}"
                        ) from exc
                    submit_next()
        except BrokenProcessPool as exc:
            shutdown_executor(executor, wait=False)
            raise RuntimeError(
                f"{stage}: a worker process exited abruptly while processing "
                f"{len(selected)} bricks. This often means the OS killed a worker "
                "for OOM; reduce --processes."
            ) from exc
        except Exception:
            shutdown_executor(executor, wait=False)
            raise
        else:
            shutdown_executor(executor, wait=True)

    elapsed = time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start))
    print(f"{stage} Done! {elapsed} {counts}")


def run_custom_stage(stage: str, args: argparse.Namespace, bricks: Table, path: str) -> None:
    star_ra, star_dec, star_radius, rectangles = parse_custom_mask(path)
    print(
        f"{stage}: circular masks={len(star_ra)} rectangular masks={len(rectangles)}"
    )
    process_bricks(
        stage=stage,
        args=args,
        bricks=unique_bricks(bricks),
        worker=worker_component,
        state={
            "star_ra": star_ra,
            "star_dec": star_dec,
            "star_radius": star_radius,
            "star_sky": skycoord_from_arrays(star_ra, star_dec),
            "rectangles": rectangles,
        },
    )


def run_star_stage(
    stage: str,
    args: argparse.Namespace,
    bricks: Table,
    catalog_path: str,
    radius_field_template: str | None,
) -> None:
    if radius_field_template is None:
        radius_column = "radius"
    else:
        # These stages are field-specific because north/south use different radii.
        radius_column = None

    if radius_column is not None:
        star_ra, star_dec, star_radius = read_star_catalog(catalog_path, radius_column)
        print(f"{stage}: stars with positive radius={len(star_ra)}")
        process_bricks(
            stage=stage,
            args=args,
            bricks=bricks,
            worker=worker_component,
            state={
                "star_ra": star_ra,
                "star_dec": star_dec,
                "star_radius": star_radius,
                "star_sky": skycoord_from_arrays(star_ra, star_dec),
                "rectangles": np.empty((0, 4), dtype="f8"),
            },
        )
        return

    for field in args.fields:
        field_bricks = bricks[np.asarray(bricks["field"]).astype(str) == field]
        if len(field_bricks) == 0:
            continue
        radius_column = radius_field_template.format(field=field)
        star_ra, star_dec, star_radius = read_star_catalog(catalog_path, radius_column)
        print(f"{stage} {field}: stars with positive radius={len(star_ra)}")
        process_bricks(
            stage=stage,
            args=args,
            bricks=field_bricks,
            worker=worker_component,
            state={
                "star_ra": star_ra,
                "star_dec": star_dec,
                "star_radius": star_radius,
                "star_sky": skycoord_from_arrays(star_ra, star_dec),
                "rectangles": np.empty((0, 4), dtype="f8"),
            },
        )


def run_combine_stage(stage: str, args: argparse.Namespace, bricks: Table) -> None:
    process_bricks(
        stage=stage,
        args=args,
        bricks=bricks,
        worker=worker_combine,
        state={},
    )


def run_stage(stage: str, args: argparse.Namespace, bricks: Table) -> None:
    print(f"\n=== {stage} ===")
    if stage == "all_custom":
        run_custom_stage(stage, args, bricks, args.all_custom_mask_file)
    elif stage == "elg_custom":
        run_custom_stage(stage, args, bricks, args.elg_custom_mask_file)
    elif stage == "lrg_gaia":
        run_star_stage(stage, args, bricks, args.lrg_gaia_catalog, "radius_{field}")
    elif stage == "lrg_wise":
        run_star_stage(stage, args, bricks, args.lrg_wise_catalog, None)
    elif stage == "elg_gaia":
        run_star_stage(stage, args, bricks, args.elg_gaia_catalog, "radius_{field}")
    elif stage in {"combine_lrg", "combine_elg"}:
        run_combine_stage(stage, args, bricks)
    else:
        raise ValueError(f"Unknown stage {stage}")


def print_configuration(args: argparse.Namespace) -> None:
    print("Input catalogs:")
    print(f"  all_custom_mask_file: {args.all_custom_mask_file}")
    print(f"  elg_custom_mask_file: {args.elg_custom_mask_file}")
    print(f"  lrg_gaia_catalog: {args.lrg_gaia_catalog}")
    print(f"  elg_gaia_catalog: {args.elg_gaia_catalog}")
    print(f"  lrg_wise_catalog: {args.lrg_wise_catalog}")
    print("Output directories:")
    print(f"  all_custom_dir: {args.all_custom_dir}")
    print(f"  elg_custom_dir: {args.elg_custom_dir}")
    print(f"  lrg_gaia_dir: {args.lrg_gaia_dir}")
    print(f"  lrg_wise_dir: {args.lrg_wise_dir}")
    print(f"  elg_gaia_dir: {args.elg_gaia_dir}")
    print(f"  lrg_output_dir: {args.lrg_output_dir}")
    print(f"  elg_output_dir: {args.elg_output_dir}")


def main() -> None:
    args = parse_args()
    stages = STAGE_ORDER if "all" in args.stages else args.stages
    print_configuration(args)
    bricks = read_bricks(args)
    bricks = select_debug_bricks(bricks, args)
    for stage in stages:
        run_stage(stage, args, bricks)


if __name__ == "__main__":
    main()
