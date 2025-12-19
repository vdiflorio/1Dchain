#!/bin/bash

#SBATCH -o output/output_%x_%j.txt
#SBATCH -e output/errors_%x_%j.txt

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=cpu_sapphire
#SBATCH --mem=1G
#SBATCH --time=00:20:59  

N_VALUES=(50 100 130 170 230 270 350)

for N in "${N_VALUES[@]}"; do
    echo "Submitting job for N=$N"
    sbatch init_cond.sh $N
done