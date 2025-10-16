# I'd like to make a run script to run main.exe with mpiexec to easily specify number of processors as an input argument from command line.
# Here's a simple example of such a script:
#!/bin/bash
# Usage: ./run.sh [num_procs]
# Default to 4 processors if not specified
if [ -n "$1" ]; then
  NUM_PROCS=$1
else
  NUM_PROCS=5
fi  
mpiexec -n $NUM_PROCS ./main.exe