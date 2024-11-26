#!/bin/bash

#SBATCH -o output/output_%x_pid_%j.txt
#SBATCH -e output/errors_%x_pid_%j.txt

#SBATCH --nodes=4 
#SBATCH --ntasks=160
#SBATCH --cpus-per-task=1 
#SBATCH --partition=long_cpu 
#SBATCH --mem=64G
#SBATCH --time=16:00:00  

module use /opt/share/sw2/modules/all
module load OpenMPI/4.1.6-GCC-13.2.0

make clean all

# Loop over argument combinations
# for N in 30 50 70; do         # Adjust these values as needed
    for Tr in 2.0 2.5 3.0; do # Adjust these values as needed
      echo "Running with N=30, Tl=1, Tr=$Tr"
      mpirun -np 2 ./fput 30 1 $Tr
    done
# done