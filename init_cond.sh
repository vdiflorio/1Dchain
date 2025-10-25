#!/bin/bash

#SBATCH -o output/output_%x_pid_%j.txt
#SBATCH -e output/errors_%x_pid_%j.txt

#SBATCH --nodes=3
#SBATCH --ntasks=190
#SBATCH --partition=cpu_sapphire
#SBATCH --mem=90G
#SBATCH --time=23:59:59  

module load /share/apps/legion-modulefiles/openmpi/5.0.7_gcc12 


JSON_FILE="parametri_simu.json"

for i in {0..9}; do
  echo "Esecuzione numero: $i"
for N in  250 500; do #100 200 300 400 500 750 1000; do 
for grad in 0.0001; do
    
    
    python3 scriptino.py "$JSON_FILE" "$N" "$grad" "$i"
        # Controlla se la modifica è riuscita
        if [ $? -ne 0 ]; then
          echo "Errore durante la modifica del file JSON. Esco."
          exit 1
        fi
        echo "File JSON modificato con successo."
        mpirun --bind-to none -np 190 ./fput  $JSON_FILE

done
done
done
