#!/bin/bash


# Percorso al file JSON e al programma
JSON_FILE="parametri.json"

for N in $(seq 30 10 30); do 
    for grad in $(awk 'BEGIN{for(i=0.01;i<=0.05;i+=0.01)printf "%.2f ", i}'); do
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
    tar -czvf single_data/compressed_mil_N_${N}.tar.gz single_data/ttcf_mil_N_${N}_Tr_*.dat
    rm single_data/ttcf_mil_N_${N}_Tr_*.dat
done