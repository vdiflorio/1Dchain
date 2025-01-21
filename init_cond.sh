#!/bin/bash

#SBATCH -o output/output_%x_pid_%j.txt
#SBATCH -e output/errors_%x_pid_%j.txt

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --partition=disma
#SBATCH --mem=32G
#SBATCH --time=1-23:00:00  



JSON_FILE="parametri.json"


for N in 50; do #100 200 300 400 500 750 1000; do 
#for grad in 0.0001 0.001 0.01 0.1 1; do
    
    echo "Running with N=$N, Tr=$Tr"
    python3 scriptino.py "$JSON_FILE" "$N" "0"
        # Controlla se la modifica è riuscita
        if [ $? -ne 0 ]; then
          echo "Errore durante la modifica del file JSON. Esco."
          exit 1
        fi
        mpirun -np 1 ./fput  $JSON_FILE

#done
done

