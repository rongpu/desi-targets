#!/usr/bin/env bash
#SBATCH --qos=regular
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --constraint=cpu
#SBATCH --time=12:00:00
#SBATCH --account=desi

set -euo pipefail

if [[ -n "${DR11MASK_DIR:-}" ]]; then
  cd "$DR11MASK_DIR"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/payload.sh" ]]; then
  cd "$SLURM_SUBMIT_DIR"
else
  cd "$(dirname "$0")"
fi

export PROCESSES=${PROCESSES:-128}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}
NNODES=${SLURM_NNODES:-4}

run_stage() {
  local tasks_file=$1
  local stage=$2
  echo "=== $stage ==="
  srun --wait=0 \
    --nodes="$NNODES" \
    --ntasks="$NNODES" \
    --ntasks-per-node=1 \
    --cpu-bind=cores \
    ./payload.sh "$tasks_file" "$stage"
}

run_stage custom_tasks.txt all_custom
run_stage custom_tasks.txt elg_custom
run_stage tasks.txt lrg_gaia
run_stage tasks.txt lrg_wise
run_stage tasks.txt elg_gaia
run_stage tasks.txt combine_lrg
run_stage tasks.txt combine_elg
