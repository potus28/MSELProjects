#!/bin/bash
#SBATCH --job-name=@@@
#SBATCH --nodes=5
#SBATCH --ntasks=200
#SBATCH --time=48:00:00
#SBATCH --partition=400p48h
#SBATCH --qos=funded
#SBATCH --account=060803-362094
#SBATCH --mail-user=wnw36@msstate.edu
#SBATCH --mail-type=ALL

module purge
module load intel/2021.2
module load impi/2021.2

export ASE_CP2K_COMMAND="srun cp2k_shell.psmp"

rootdir=

lignin=$rootdir/Resources/xyz/lignin/GG/BO4.xyz
zeolite=$rootdir/DOE/AICHE2022/laptop/fire.traj

for i in pureSi Al1 Al4 Al4-Al1; do

    mkdir $i
    cp run.py $i
    cd $i
    ln -s ../../../../../Python/mycalculators.py mycalculators.py
    python3 run.py $i $lignin $zeolite $i
    cd ..

done
