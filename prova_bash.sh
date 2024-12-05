#!/bin/bash



#SBATCH -o output/output_%x_pid_%j.txt
#SBATCH -e output/errors_%x_pid_%j.txt

#SBATCH --nodes=4 
#SBATCH --ntasks=160
#SBATCH --cpus-per-task=1 
#SBATCH --partition=long_cpu 
#SBATCH --mem=128G
#SBATCH --time=23:00:00  

# Percorso al file JSON e al programma
JSON_FILE="parametri.json"
make clean all
for N in $(seq 50 50 210); do 
    for grad in $(awk 'BEGIN{for(i=0.01;i<=0.05;i+=0.01)printf "%.2f ", i}'); do
        echo "Running with N=$N, gradT=$grad"
        python3 scriptino.py "$JSON_FILE" "$N" "$grad"
        # Controlla se la modifica è riuscita
        if [ $? -ne 0 ]; then
          echo "Errore durante la modifica del file JSON. Esco."
          exit 1
        fi
        mpirun -np 160 ./fput  $JSON_FILE
    done
    # Compress all files for the current N after completing the inner loop
    tar -czvf script/compressed_archive_N_${N}.tar.gz script/ttcf_mil_N_${N}_Tr_*.dat
    rm script/ttcf_mil_N_${N}_Tr_*.dat
    python3 costante.py ${N}
done