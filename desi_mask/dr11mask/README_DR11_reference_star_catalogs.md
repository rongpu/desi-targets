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


## Per-brick pixel masks

After the reference catalogs exist, build the DR11 per-brick pixel masks with:

```bash
python build_dr11_pixel_masks.py --help
```

The pixel-mask builder keeps the DR9 mask definitions and directory layout, but
uses Cartesian unit-vector WCS interpolation instead of scalar RA/Dec
interpolation. This avoids the RA-wrap and celestial-pole failure mode in the
DR9 rasterization scripts.

Typical command for one field/task split on a Perlmutter CPU node:

```bash
python build_dr11_pixel_masks.py \
  --release-root /dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr11 \
  --reference-dir /global/cfs/cdirs/desi/users/$USER/desi_mask/dr11_reference_stars \
  --output-root /global/cfs/cdirs/desi/users/$USER/desi_mask/dr11_pixel_masks \
  --fields south \
  --n-task 4 \
  --task-id 0 \
  --processes 128 \
  --stages lrg_gaia
```

The stage order is:

```text
all_custom
elg_custom
lrg_gaia
lrg_wise
elg_gaia
combine_lrg
combine_elg
```

For DR9-style multi-node dispatch with GNU parallel, submit:

```bash
sbatch launch_dr11_pixel_masks.sh
```

`launch_dr11_pixel_masks.sh` requests four Perlmutter CPU nodes by default.
It runs one `srun` step per stage, so custom/component stages finish
across all nodes before the combine stages start. For each stage, `payload.sh`
runs once per node, loads GNU parallel, assigns rows from `custom_tasks.txt`
(custom stages) or `tasks.txt` (field-specific component/combine stages) by
`SLURM_NODEID`, and calls `run_dr11_pixel_masks.sh`. Each split then uses
Python multiprocessing with 128 local workers.

`run_dr11_pixel_masks.sh` can still be run directly for debugging. Set
`STAGES="lrg_gaia"` or another space-separated stage list to run selected
stages. A direct all-stage run requires `N_TASK=1`; the script refuses
`STAGES=all` with `N_TASK>1` because direct mode has no multi-node stage
barriers.

The final LRG mask bits are:

```text
0  DESI targeting MASKBITS bits 1,12,13
1  WISEM1 bits 0,1,2,3,4,6,7
2  all-tracer custom mask
3  LRG Gaia mask
4  LRG WISE mask
```

The final ELG mask bits are:

```text
0  DESI targeting MASKBITS bits 1,12,13
1  all-tracer custom mask
2  ELG Gaia mask
3  ELG-specific custom mask
```
