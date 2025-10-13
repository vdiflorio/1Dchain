#!/bin/bash

#SBATCH -o output/output_%x_pid_%j.txt
#SBATCH -e output/errors_%x_pid_%j.txt

#SBATCH --nodes=6
#SBATCH --ntasks=180
#SBATCH --partition=cpu_skylake
#SBATCH --mem=128G
#SBATCH --time=1-23:59:00  

module load /share/apps/legion-modulefiles/openmpi/5.0.7_gcc12 


JSON_FILE="parametri_test.json"


for N in 600; do #100 200 300 400 500 750 1000; do 
for grad in 0.0001; do
    
    
    python3 scriptino.py "$JSON_FILE" "$N" "$grad"
        # Controlla se la modifica è riuscita
        if [ $? -ne 0 ]; then
          echo "Errore durante la modifica del file JSON. Esco."
          exit 1
        fi
        echo "File JSON modificato con successo."
        mpirun --bind-to none -np 180 ./fput  $JSON_FILE

done
done

