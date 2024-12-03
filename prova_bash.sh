#!/bin/bash

# Controllo degli argomenti
if [ "$#" -ne 2 ]; then
  echo "Uso: $0 <Nuovo_N> <Nuovo_Tr>"
  exit 1
fi

# Valori passati come argomenti
NEW_N=$1
NEW_TR=$2

# Percorso al file JSON e al programma
JSON_FILE="parametri.json"

# Modifica il JSON usando lo script Python
python3 modify_n_tr.py "$JSON_FILE" "$NEW_N" "$NEW_TR"

# Controlla se la modifica è riuscita
if [ $? -ne 0 ]; then
  echo "Errore durante la modifica del file JSON. Esco."
  exit 1
fi


for N in $(seq 30 10 40); do 
    for grad in $(seq 1.001 0.002 1.1); do
        echo "Running with N=$N, gradT=$grad"
        python scriptino.py "$JSON_FILE" "$N" "$grad"
        mpirun -np 4 ./fput  $JSON_FILE
    done
    # Compress all files for the current N after completing the inner loop
    tar --remove-files -czvf single_data/compressed_mil_N_${N}.tar.gz single_data/ttcf_mil_N_${N}_Tr_*.dat
done