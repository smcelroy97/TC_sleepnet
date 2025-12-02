#!/bin/bash
#SBATCH --job-name=NO_diff_test
#SBATCH -A TG-MED240050
#SBATCH -t 5:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH -o simOutput/single.run
#SBATCH -e simOutput/single.err

#SBATCH --mem=450G
#SBATCH --export=ALL
#SBATCH --partition=GPU
source ~/.bashrc
mpiexec -n 64 nrniv -mpi -python bazh_net.py