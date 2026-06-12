# DR11 per-brick pixel masks

This document describes the DR11 per-brick pixel-mask builder. The reference
catalog inputs are built separately; see `README_DR11_reference_star_catalogs.md`.

The main entry point is:

```bash
python build_dr11_pixel_masks.py --help
```

The pixel-mask builder keeps the DR9 mask definitions and directory layout, but
uses Cartesian unit-vector WCS interpolation instead of scalar RA/Dec
interpolation. This avoids the RA-wrap and celestial-pole failure mode in the
DR9 rasterization scripts.

## Typical command

Typical command for one field/task split on a Perlmutter CPU node:

```bash
python build_dr11_pixel_masks.py \
  --release-root /dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr11 \
  --reference-dir /global/cfs/cdirs/desi/users/$USER/desi_mask/dr11_reference_stars \
  --output-root /global/cfs/cdirs/desi/users/$USER/desi_mask/dr11_pixel_masks \
  --fields south \
  --n-task 4 \
  --task-id 0 \
  --processes 32 \
  --stages lrg_gaia
```

## Stage order

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

## Multi-node dispatch

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
Python multiprocessing with 32 local workers.

`run_dr11_pixel_masks.sh` can still be run directly for debugging. Set
`STAGES="lrg_gaia"` or another space-separated stage list to run selected
stages. A direct all-stage run requires `N_TASK=1`; the script refuses
`STAGES=all` with `N_TASK>1` because direct mode has no multi-node stage
barriers.

## Output mask bits

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
