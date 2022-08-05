#!/bin/bash

for t in 300.0; do
	mkdir -p $t
	for mu in 50.0 45.0 40.0 35.0 30.0 25.0 20.0 15.0 10.0 5.0 1.0; do
		mkdir -p $t/$mu
			cp -r init_files/* $t/$mu
			cd $t/$mu
			sed -i "s/temptemptemp/$t/g" gcmc.eq.inp
			sed -i "s/mumumu/$mu/g" gcmc.eq.inp
			$cassandra gcmc.eq.inp >> cassandra.out
			cd ../..
	done
done
