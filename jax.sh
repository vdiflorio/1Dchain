#!/bin/bash

#SBATCH -o output/output_%x_pid_%j.txt
#SBATCH -e output/errors_%x_pid_%j.txt

#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --time=23:59:00  
#SBATCH --partition=gpu_a40
#SBATCH --gres=gpu:1 

conda deactivate

cd training
conda activate jax-env

python3 test_jax.py