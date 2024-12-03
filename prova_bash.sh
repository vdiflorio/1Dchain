#!/bin/bash


# Percorso al file JSON e al programma
JSON_FILE="parametri.json"

for N in $(seq 30 10 40); do 
    for grad in $(seq 0.01 0.05 0.01); do
        echo "Running with N=$N, gradT=$grad"
        python scriptino.py "$JSON_FILE" "$N" "$grad"
        # Controlla se la modifica è riuscita
        if [ $? -ne 0 ]; then
          echo "Errore durante la modifica del file JSON. Esco."
          exit 1
        fi
        mpirun -np 4 ./fput  $JSON_FILE
    done
    # Compress all files for the current N after completing the inner loop
    tar --remove-files -czvf single_data/compressed_mil_N_${N}.tar.gz single_data/ttcf_mil_N_${N}_Tr_*.dat
done