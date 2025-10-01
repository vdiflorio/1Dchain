#!/bin/bash

#SBATCH -o output/output_%x_pid_%j.txt
#SBATCH -e output/errors_%x_pid_%j.txt

#SBATCH --nodes=4
#SBATCH --ntasks=128
#SBATCH --partition=cpu_sapphire
#SBATCH --mem=32G
#SBATCH --time=23:00:00  



JSON_FILE="parametri_simu.json"


for N in 30; do #100 200 300 400 500 750 1000; do 
for grad in 0.0001; do
    
    
    python3 scriptino.py "$JSON_FILE" "$N" "$grad"
        # Controlla se la modifica è riuscita
        if [ $? -ne 0 ]; then
          echo "Errore durante la modifica del file JSON. Esco."
          exit 1
        fi
        echo "File JSON modificato con successo."
        mpirun -np 128 ./fput  $JSON_FILE

done
done

