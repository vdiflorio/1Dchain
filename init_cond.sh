#!/bin/bash

#SBATCH -o output/output_%x_pid_%j.txt
#SBATCH -e output/errors_%x_pid_%j.txt

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --partition=disma 
#SBATCH --mem=128G
#SBATCH --time=1-23:00:00  



JSON_FILE="parametri.json"

for N in 50 110 250 500 1000; do

        grad=$(echo "scale=8; 127 / $N" | bc)
        echo "Running with N=$N, gradT=$grad"
        python3 scriptino.py "$JSON_FILE" "$N" "$grad"
        # Controlla se la modifica è riuscita
        if [ $? -ne 0 ]; then
          echo "Errore durante la modifica del file JSON. Esco."
          exit 1
        fi
        mpirun -np 1 ./fput  $JSON_FILE


done

