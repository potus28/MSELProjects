#!/bin/bash
export OMP_NUM_THREADS=2
mkdir -p simulation
for t in 300.0; do

	mkdir -p simulation/$t
	for mu in 50.0 45.0 40.0 35.0 30.0 25.0 20.0 15.0; do

            mkdir -p simulation/$t/$mu
			cp init_files/gcmc.prd.inp simulation/$t/$mu
			cd simulation/$t/$mu
			sed -i "s/temptemptemp/$t/g" gcmc.prd.inp
			sed -i "s/mumumu/$mu/g" gcmc.prd.inp

            sed -i '/nmols/a mass_density' gcmc.prd.inp
            sed -i '/nmols/a density' gcmc.prd.inp
            sed -i '/nmols/a volume' gcmc.prd.inp
            sed -i '/nmols/a pressure' gcmc.prd.inp
            sed -i '/nmols/a enthalpy' gcmc.prd.inp

            cassandra gcmc.prd.inp >> cassandra.out &
			cd ../../..
	done
done

