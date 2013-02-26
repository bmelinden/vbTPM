#!/bin/bash

# Unix/linux shell script to start repeated VB7 computations in the
# background. If the runinput parameter one_at_a_time=true, then this
# means Matlab is restarted between simulation runs, which is a
# workaround for Matlab's bad memory management. 

# How to use:

# 1) Make a copy of this file, and adjust the name of the runinput
# file (myruninputfile.m) and the number of iterations (second number
# in the expression (1..50). The copy should be in the same directory
# as the runinput file.

# 2) Add the VB7 directories in the Matlab startup path, e.g., by
# editing startup.m (see Matlab documentation).

# 2) Make the new copy executable, by issuing something like
# > chmod +x myrun.sh

# 3) To run jobs over a remote ssh connection, start jobs in the
# background using nohup to avoid termination at logout, and pipe all
# output to a log file.

# In a bash shell, issue something like
# nohup ./myrun.sh &> log1 & 
# nohup ./myrun.sh &> log2 & 
# nohup ./myrun.sh &> log3 & 
# etc. Do not start more jobs than you have cores on the machine.

# These runs should continue to run even if you log out, and you can
# monitor the progress by inspecting the log files. 

# the actual script:
for i in {1..50}
do
    nohup echo "VB7_batch_run('myruninputfile.m'),exit" | nice -n 15 matlab -nosplash -nodesktop
done

