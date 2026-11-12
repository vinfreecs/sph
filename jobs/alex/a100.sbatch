#!/bin/bash -l
# NHR@FAU Alex, one NVIDIA A100. Build first (see README), then from the repo root:
#   sbatch jobs/alex/a100.sbatch [sph_cuda options, e.g. --scene box --out output/a100]
#SBATCH --job-name=sph-a100
#SBATCH --gres=gpu:a100:1
#SBATCH --partition=a100
#SBATCH --time=01:00:00
#SBATCH --export=NONE

unset SLURM_EXPORT_ENV
module load cuda

cd "$SLURM_SUBMIT_DIR"
./build/sph_cuda "$@"
