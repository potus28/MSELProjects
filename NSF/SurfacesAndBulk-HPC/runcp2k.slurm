#!/bin/bash
#SBATCH --job-name=MYCP2K
#SBATCH --output=MYCP2K.o%j
#SBATCH --error=MYCP2K.e%j
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --partition=wholenode
#SBATCH --time=30:00:00
#SBATCH --account=che140141
ulimit -s unlimited
export OMP_NUM_THREADS=1
echo "SLURM_NTASKS: " $SLURM_NTASKS

export ASE_CP2K_COMMAND="mpirun -np $SLURM_NTASKS cp2k_shell.psmp"

metals=("Alpha-Mo2C_mp-1552" "Beta-Mo2C_mp-1221498")
functionals=("optB88-vdw" "optB86B-vdw" "optPBE-vdw" "rVV10" "PBE+D3")
miller=("001" "101" "121"  "120"  "002")

for m in ${metals[@]}; do
	for f in ${functionals[@]}; do
		python3 run-bulk.py $m $f
		for idx in ${miller[@]}; do
			python3 run-surface.py $m $f $idx
		done
	done
done

