#!/bin/bash -l
# NHR@FAU Alex, one NVIDIA A40. Build first (see README), then from the repo root:
#   sbatch jobs/alex/a40.sbatch [sph_cuda options, e.g. --scene box --out output/a40]
#SBATCH --job-name=sph-a40
#SBATCH --gres=gpu:a40:1
#SBATCH --partition=a40
#SBATCH --time=01:00:00
#SBATCH --export=NONE

unset SLURM_EXPORT_ENV
module load cuda

cd "$SLURM_SUBMIT_DIR"
./build/sph_cuda "$@"
