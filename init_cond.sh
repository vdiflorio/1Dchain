#!/bin/bash

#SBATCH -o output/output_%x_pid_%j.txt
#SBATCH -e output/errors_%x_pid_%j.txt

#SBATCH --nodes=3
#SBATCH --ntasks=190
#SBATCH --partition=cpu_sapphire
#SBATCH --mem=90G
#SBATCH --time=23:59:59  

module load /share/apps/legion-modulefiles/openmpi/5.0.7_gcc12 

TEMPLATE_JSON="parametri_simu.json"

for N in 250; do  # 100 200 300 400 500 750 1000
for i in {0..15}; do
  echo "Esecuzione numero: $i"

  for grad in 0.0001; do

    # Crea una copia temporanea del file JSON
    JSON_COPY="parametri_simu_${N}_${grad}_${i}.json"
    cp "$TEMPLATE_JSON" "$JSON_COPY"

    echo "Utilizzo il file JSON copiato: $JSON_COPY"

    # Modifica la copia tramite lo script Python
    python3 scriptino.py "$JSON_COPY" "$N" "$grad" "$i"
    if [ $? -ne 0 ]; then
      echo "Errore durante la modifica del file JSON. Esco."
      exit 1
    fi
    echo "File JSON modificato con successo."

    # Esegui il programma MPI usando la copia
    mpirun --bind-to none -np 190 ./fput "$JSON_COPY"

    # (Opzionale) Rimuovi la copia per non lasciare file temporanei
    rm -f "$JSON_COPY"

  done
done
done
