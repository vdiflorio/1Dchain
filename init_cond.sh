#!/bin/bash

#SBATCH -o output/output_%x_pid_%j.txt
#SBATCH -e output/errors_%x_pid_%j.txt

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --partition=long_cpu
#SBATCH --mem=128G
#SBATCH --time=1-23:00:00  



JSON_FILE="parametri_init.json"

Tr=128
for N in 50 100 150 500 1000; do 
    #for grad in $(awk 'BEGIN{for(i=0.001;i<=0.02;i+=0.002)printf "%.4f ", i}'); do
    
    echo "Running with N=$N, Tr=$Tr"
    python3 scriptino.py "$JSON_FILE" "$N" "$Tr"
        # Controlla se la modifica è riuscita
        if [ $? -ne 0 ]; then
          echo "Errore durante la modifica del file JSON. Esco."
          exit 1
        fi
        mpirun -np 1 ./fput  $JSON_FILE


done

