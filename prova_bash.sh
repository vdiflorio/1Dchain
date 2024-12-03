#!/bin/bash

# Percorso al file JSON e al programma
JSON_FILE="config.json"
PROGRAM="./programma_eseguibile"

# Modifica il JSON usando lo script Python
python3 scriptino.py "$JSON_FILE"

# Controlla se la modifica è andata a buon fine
if [ $? -ne 0 ]; then
  echo "Errore durante la modifica del file JSON. Esco."
  exit 1
fi

# Esegui il programma che utilizza il file JSON
$PROGRAM "$JSON_FILE"
