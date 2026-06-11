# DR11 reference star catalogs

This directory contains a DR11-oriented replacement for the DR9 reference-star
catalog scripts in `desi_mask/reference/`, `desi_mask/lrg_mask/`, and
`desi_mask/elg_mask/`.

The main entry point is:

```bash
python build_dr11_reference_star_catalogs.py --help
```

The script keeps the DR9 mask-radius prescriptions, but makes the footprint and
output naming configurable. By default it assumes the Legacy Survey reduction is
available at:

```text
/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr11
```

and that the brick files are:

```text
{release_root}/{field}/survey-bricks-dr11-{field}.fits.gz
```

where `field` is `north` or `south`.

## Typical command

```bash
python build_dr11_reference_star_catalogs.py \
  --release-root /dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr11 \
  --output-dir /global/cfs/cdirs/desi/users/$USER/desi_mask/dr11_reference_stars \
  --stages all
```

## Important choices

- The Gaia source input remains Gaia EDR3 healpix files by default.
- Gaia healpix trimming reads are parallelized by default with
  `--gaia-read-processes` capped at 8 workers; pass `--gaia-read-processes 1`
  to force serial reads.
- The LRG and ELG Gaia radius-vs-magnitude relations are copied from the DR9
  scripts.
- The LRG WISE radius-vs-magnitude relation is copied from the DR9 script.
- The old DR9 Gaia supplement is optional. For new DR11-only sky, there is no
  DR9 supplement by construction. Pass `--old-gaia-mask PATH` only if you want
  to preserve the old-overlap supplement behavior.

## Outputs

The default output names are:

```text
tycho2-reference-dr11.fits
gaia_edr3_g_18_dr11.fits
gaia_edr3_g_14_pm_dr11.fits
gaia_edr3_g_18_predict_decam_dr11.fits
gaia_reference_dr11.fits
gaia_reference_suppl_dr11.fits        # only when --old-gaia-mask is supplied
gaia_lrg_mask_dr11_v1.fits
gaia_elg_mask_dr11_v1.fits
w1_bright-2mass-dr11.fits
w1_bright-13.3-dr11.fits
w1_bright-2mass-13.3-dr11.fits
w1_bright-2mass-lrg_mask_dr11_v1.fits
```

