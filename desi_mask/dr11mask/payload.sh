#!/usr/bin/env bash
set -euo pipefail

tasks_file=${1:-tasks.txt}
stage=${2:-all}

if [[ ! -f "$tasks_file" ]]; then
  echo "tasks file not found: $tasks_file" >&2
  exit 1
fi
if [[ -z "${SLURM_NODEID:-}" ]]; then
  echo "need SLURM_NODEID set" >&2
  exit 1
fi
if [[ -z "${SLURM_NNODES:-}" ]]; then
  echo "need SLURM_NNODES set" >&2
  exit 1
fi

if type module >/dev/null 2>&1; then
  module load parallel
fi

if [[ -n "${DESI_ENV_VERSION:-}" ]]; then
  source /dvs_ro/cfs/cdirs/desi/software/desi_environment.sh "$DESI_ENV_VERSION"
fi

if ! command -v parallel >/dev/null 2>&1; then
  echo "GNU parallel is not available; load it or run module load parallel" >&2
  exit 1
fi

export PROCESSES=${PROCESSES:-128}
export PYTHON=${PYTHON:-python}
export STAGES=${stage}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}

awk -v NNODE="$SLURM_NNODES" -v NODEID="$SLURM_NODEID" '
  NF && $1 !~ /^#/ {
    task += 1
    if (((task - 1) % NNODE) == NODEID) print
  }
' "$tasks_file" | \
parallel --jobs 1 --colsep ' ' \
  'FIELD={1} N_TASK={2} TASK_ID={3} bash run_dr11_pixel_masks.sh'
