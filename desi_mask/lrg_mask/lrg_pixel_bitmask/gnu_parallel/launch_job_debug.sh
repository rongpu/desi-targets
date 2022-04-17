#!/bin/bash
#SBATCH --qos=debug
#SBATCH --nodes=3
#SBATCH --ntasks-per-node 1
#SBATCH --constraint=haswell
#SBATCH --time=00:30:00
#SBATCH --account=desi

srun --wait=0 payload.sh tasks.txt
