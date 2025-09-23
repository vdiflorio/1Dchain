#!/bin/bash


#SBATCH -A dcarbone
#SBATCH -o output/output_%x_pid_%j.txt
#SBATCH -e output/errors_%x_pid_%j.txt

#SBATCH --ntasks=100
#SBATCH --cpus-per-task=1 
#SBATCH --mem=64G
#SBATCH --time=23:59:00

# Percorso al file JSON e al programma
JSON_FILE="parametri_simu.json"
#make clean all


for N in 500 600 750; do 
for grad in 0.001 0.0001; do
    echo "Running with N=$N, grad=$grad"
    
    SECONDS=0  # Reset timer
    
    python3 scriptino.py "$JSON_FILE" "$N" "$grad"
    srun ./fput
    
    duration=$SECONDS
    echo "Loop completed in $((duration / 60)) min and $((duration % 60)) sec"
done
done
