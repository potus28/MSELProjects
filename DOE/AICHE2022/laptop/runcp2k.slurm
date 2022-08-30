#!/bin/bash
#SBATCH --job-name=mycp2k
#SBATCH --output=mycp2k.o%j      
#SBATCH --error=mycp2k.e%j 
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=40
#SBATCH --time=04:00:00
#SBATCH --partition=400p48h
#SBATCH --qos=funded
#SBATCH --account=060803-362129
#SBATCH --mail-user=wnw36@msstate.edu
#SBATCH --mail-type=ALL
#SBATCH --exclusive
ulimit -s unlimited
export OMP_NUM_THREADS=1

exe=/work/hpc/users/wnw36/Programs/cp2k/2022.1/exe/Linux-gnu-x86_64/cp2k_shell.psmp

export ASE_CP2K_COMMAND="mpirun -np $SLURM_NTASKS $exe"

module purge
module load gcc/11.3.0
module load openmpi/4.0.4


python3 AICHE.py
