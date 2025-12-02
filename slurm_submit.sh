#!/bin/bash

#source ~/.bashrc

#module purge
#module load shared
#module load cpu/0.17.3b  gcc/10.2.0/npcyll4 openmpi/4.1.1/ygduf2r
#module load sdsc
#module load cpu
#conda avtivate py310

#SBATCH --job-name=NO_diff_test
#SBATCH -A TG-MED240050
#SBATCH -t 3:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH -o simOutput/single.run
#SBATCH -e simOutput/single.err

#SBATCH --mem=240G
#SBATCH --export=ALL
#SBATCH --partition=compute

#cd /home/smcelroy/TC_sleepnet

mpirun -n 64 nrniv -python -mpi bazh_net.py