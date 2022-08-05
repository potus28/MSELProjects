#!/bin/bash

for t in 300.0; do
	mkdir -p $t
	for mu in 50.0 45.0 40.0 35.0 30.0 25.0 20.0 15.0; do
		mkdir -p $t/$mu
			cp init_files/gcmc.prd.inp $t/$mu
			cd $t/$mu
			sed -i "s/temptemptemp/$t/g" gcmc.prd.inp
			sed -i "s/mumumu/$mu/g" gcmc.prd.inp
			$cassandra gcmc.prd.inp >> cassandra.out
			cd ../..
	done
done
