#!/bin/bash
#SBATCH --job-name=MYCP2K
#SBATCH --output=MYCP2K.o%j
#SBATCH --error=MYCP2K.e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --partition=wholenode
#SBATCH --time=00:30:00
#SBATCH --account=che140141

ulimit -s unlimited
export OMP_NUM_THREADS=1
echo "SLURM_NTASKS: " $SLURM_NTASKS

export ASE_CP2K_COMMAND="mpirun -np $SLURM_NTASKS cp2k_shell.psmp"

libfile=/home/woodrowwilson/Documents/Projects/MSELProjects/Python/mycalculators.py
metals=(Alpha-Mo2C_mp-1552 Beta-Mo2C_mp-1221498)
functionals=(rVV10)
probe=(FAL.xyz)

m=Alpha-Mo2C_mp-1552
p=FAL
f=rVV10

mkdir -p $f/$m/$p
cp MD.py $f/$m/$p
cd $f/$m/$p

ln -s $libfile mycalculators.py
python3 MD.py ${p}.xyz $m $f

