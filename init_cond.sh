#!/bin/bash

#SBATCH -o output/output_%x_pid_%j.txt
#SBATCH -e output/errors_%x_pid_%j.txt

#SBATCH --nodes=4 
#SBATCH --ntasks=160
#SBATCH --cpus-per-task=1 
#SBATCH --partition=long_cpu 
#SBATCH --mem=128G
#SBATCH --time=16:00:00  

module use /opt/share/sw2/modules/all
module load OpenMPI/4.1.6-GCC-13.2.0

make clean all



mpirun -np 1 ./fput 80 1 1 y
