#!/bin/bash

#SBATCH -o output/output_%x_pid_%j.txt
#SBATCH -e output/errors_%x_pid_%j.txt

#SBATCH --nodes=4 
#SBATCH --ntasks=160
#SBATCH --cpus-per-task=1 
#SBATCH --partition=long_cpu 
#SBATCH --mem=128G
#SBATCH --time=23:00:00  

module use /opt/share/sw2/modules/all
module load OpenMPI/4.1.6-GCC-13.2.0

make clean all

for N in $(seq 40 10 200); do 
    for Tr in $(seq 1.001 0.002 1.1); do
        echo "Running with N=$N, Tl=1, Tr=$Tr"
        mpirun -np 160 ./fput $N 1 $Tr n
    done

    # Compress all files for the current N after completing the inner loop
    tar --remove-files -czvf single_data/compressed_mil_N_${N}.tar.gz single_data/ttcf_mil_N_${N}_Tr_*.dat
done