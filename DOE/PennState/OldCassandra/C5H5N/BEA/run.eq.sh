#!/bin/bash

export OMP_NUM_THREADS=2

mkdir -p simulation
for t in 300.0; do
	mkdir -p simulation/$t
	for mu in 50.0 45.0 40.0 35.0 30.0 25.0 20.0 15.0; do
		mkdir -p simulation/$t/$mu
			cp -r init_files/* simulation/$t/$mu
			cd simulation/$t/$mu
			sed -i "s/temptemptemp/$t/g" gcmc.eq.inp
			sed -i "s/mumumu/$mu/g" gcmc.eq.inp

            sed -i '/nmols/a mass_density' gcmc.eq.inp
            sed -i '/nmols/a density' gcmc.eq.inp
            sed -i '/nmols/a volume' gcmc.eq.inp
            sed -i '/nmols/a pressure' gcmc.eq.inp
            sed -i '/nmols/a enthalpy' gcmc.eq.inp


			cassandra.exe gcmc.eq.inp >> cassandra.out &
            cd ../../..
	done
done

wait
