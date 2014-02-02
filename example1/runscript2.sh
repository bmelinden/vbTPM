#!/bin/bash

# This is a bash-script used to run analysis on Ubuntu linux
# systems. VB7_batch_run.m keeps tracks of which trajectories are
# done, so simple parallellization is achieved by running several
# scripts simultaneously (one for each CPU). On other systems, you
# will have to create your own script.

for RI in runinput2.m # loop over multiple runinput files (ony one here).
do
    for i in {1..5}   # number of trajectories to analyze at each call
		      # to this script
    do
	echo "########## file" $RI", round" $i 
	# main call to matlab to perform one round of analysis
	echo "VB7_batch_run('${RI}'),exit" | nice -n 15 matlab -nosplash -nodesktop
	# the reason for this construction (one analysis round, then
	# restart matlab) is our experience that matlab hoards memory
	# when large data sets are loaded and discarded several times,
	# and restarting matlab is the simplest way to release that
	# memory.
    done
done

echo checking out from $(hostname)
