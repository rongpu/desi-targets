#!/usr/bin/env bash
set -euo pipefail

# cd /home/jon/temp/codex_play/dr11mask

export RELEASE=dr11
export RELEASE_ROOT=/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr11
export OUT=/global/cfs/cdirs/desi/users/${USER}/desi_mask/dr11_reference_stars

COMMON=(
  --release "$RELEASE"
  --release-root "$RELEASE_ROOT"
  --output-dir "$OUT"
)

python build_dr11_reference_star_catalogs.py "${COMMON[@]}" --stages trim_tycho2
python build_dr11_reference_star_catalogs.py "${COMMON[@]}" --stages trim_gaia_g18
python build_dr11_reference_star_catalogs.py "${COMMON[@]}" --stages trim_gaia_g14_pm
python build_dr11_reference_star_catalogs.py "${COMMON[@]}" --stages predict_decam
python build_dr11_reference_star_catalogs.py "${COMMON[@]}" --stages build_gaia_reference

# Optional old-DR9 supplement for old-overlap compatibility.
# Uncomment this block if you want to include it.
# python build_dr11_reference_star_catalogs.py "${COMMON[@]}" \
#   --old-gaia-mask /dvs_ro/cfs/cdirs/desi/users/rongpu/tmp/gaia-mask-dr9.fits \
#   --stages build_gaia_supplement

python build_dr11_reference_star_catalogs.py "${COMMON[@]}" --stages build_lrg_gaia
python build_dr11_reference_star_catalogs.py "${COMMON[@]}" --stages build_elg_gaia

python build_dr11_reference_star_catalogs.py "${COMMON[@]}" --stages trim_wise_2mass
python build_dr11_reference_star_catalogs.py "${COMMON[@]}" --stages trim_wise_faint
python build_dr11_reference_star_catalogs.py "${COMMON[@]}" --stages combine_wise
python build_dr11_reference_star_catalogs.py "${COMMON[@]}" --stages build_lrg_wise

