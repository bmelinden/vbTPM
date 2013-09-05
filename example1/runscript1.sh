#!/bin/bash
for RI in runinput1.m
do

    for i in {1..5}
    do

	echo "########## file" $RI", round" $i 
	echo "VB7_batch_run('${RI}'),exit" | nice -n 15 matlab -nosplash -nodesktop
    done
done

echo checking out from $(hostname)
