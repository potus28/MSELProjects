#!/bin/bash

metals=("Alpha-Mo2C_mp-1552" "Beta-Mo2C_mp-1221498")
functionals=("optB88-vdw" "optB86B-vdw" "optPBE-vdw" "rVV10" "PBE+D3")
miller=("001" "101" "121"  "120"  "002")

for m in ${metals[@]}; do
	for f in ${functionals[@]}; do
		echo "python3 run-bulk.py $m $f"
		for idx in ${miller[@]}; do
			echo "python3 run-surface.py $m $f $idx"
		done
	done
done

