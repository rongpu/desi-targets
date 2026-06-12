#!/usr/bin/env bash
set -euo pipefail

# Driver for one split of the DR11 per-brick pixel-mask pipeline.
# Defaults assume one Perlmutter CPU node with 128 cores.
# Run directly for one split, or let payload.sh call it from GNU parallel.

export RELEASE=${RELEASE:-dr11}
export RELEASE_ROOT=${RELEASE_ROOT:-/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr11}
export REF=/global/cfs/cdirs/desicollab/users/rongpu/desi_mask/dr11_reference_stars
export OUT=/global/cfs/cdirs/desicollab/users/rongpu/data/veto_masks/dr11

# Four-way split for a 4-node Perlmutter CPU allocation.
export FIELD=${FIELD:-south}
export N_TASK=${N_TASK:-4}
export TASK_ID=${TASK_ID:-0}
export PROCESSES=${PROCESSES:-128}
export PYTHON=${PYTHON:-python}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}

COMMON=(
  --release "$RELEASE"
  --release-root "$RELEASE_ROOT"
  --reference-dir "$REF"
  --output-root "$OUT"
  --n-task "$N_TASK"
  --task-id "$TASK_ID"
  --processes "$PROCESSES"
  # --debug-64-bricks
)
CUSTOM_COMMON=("${COMMON[@]}" --fields south north)
FIELD_COMMON=("${COMMON[@]}" --fields "$FIELD")

if [[ "${STAGES:-all}" == "all" ]]; then
  if [[ "$N_TASK" != "1" ]]; then
    echo "STAGES=all with N_TASK=$N_TASK is unsafe without stage barriers." >&2
    echo "Use launch_dr11_pixel_masks.sh, set STAGES to one stage, or set N_TASK=1 for a direct debug run." >&2
    exit 1
  fi
  STAGE_LIST=(all_custom elg_custom lrg_gaia lrg_wise elg_gaia combine_lrg combine_elg)
else
  read -r -a STAGE_LIST <<< "$STAGES"
fi

for STAGE in "${STAGE_LIST[@]}"; do
  case "$STAGE" in
    all_custom|elg_custom)
      "$PYTHON" build_dr11_pixel_masks.py "${CUSTOM_COMMON[@]}" --stages "$STAGE"
      ;;
    lrg_gaia|lrg_wise|elg_gaia|combine_lrg|combine_elg)
      "$PYTHON" build_dr11_pixel_masks.py "${FIELD_COMMON[@]}" --stages "$STAGE"
      ;;
    *)
      echo "unknown stage: $STAGE" >&2
      exit 1
      ;;
  esac
done
