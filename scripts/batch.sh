#!/bin/bash -l
#
#SBATCH --gres=gpu:a40:1
#SBATCH --time=1:00:00
#SBATCH --export=NONE

unset SLURM_EXPORT_ENV

module load cuda
cd ../build
make
./sph

