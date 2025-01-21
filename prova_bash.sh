#!/bin/bash



#SBATCH -o output/output_%x_pid_%j.txt
#SBATCH -e output/errors_%x_pid_%j.txt

#SBATCH --nodes=7
#SBATCH --ntasks=224
#SBATCH --cpus-per-task=1 
#SBATCH --partition=disma
#SBATCH --mem=32G
#SBATCH --time=1-23:00:00  

# Percorso al file JSON e al programma
JSON_FILE="parametri_simu.json"
#make clean all


for N in 600 750 800 1000; do 
for grad in 0.0001 0.001 0.01 0.1 1; do
    
    echo "Running with N=$N, Tr=$Tr"
    python3 scriptino.py "$JSON_FILE" "$N" "$grad"
        # echo "Running with N=$N, gradT=$grad"
        # python3 scriptino.py "$JSON_FILE" "$N" "$grad"
        # # Controlla se la modifica è riuscita
        # if [ $? -ne 0 ]; then
        #   echo "Errore durante la modifica del file JSON. Esco."
        #   exit 1
        # fi
        mpirun -np 224 ./fput  $JSON_FILE
    #done
    # Compress all files for the current N after completing the inner loop
    # tar -czvf short/compressed_archive_N_${N}.tar.gz short/ttcf_mil_N_${N}_Tr_*.dat
    # rm short/ttcf_mil_N_${N}_Tr_*.dat
    # python3 costante.py ${N}
done
done