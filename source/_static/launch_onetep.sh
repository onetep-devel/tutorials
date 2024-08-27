#!/bin/bash

############################################################
# This simple script is what ASE will use to start ONETEP. #
# You will need to adjust it to your environment.          #
############################################################

# Select the number of processes and threads.
# ONETEP will run on nprocesses * nthreads CPU cores.
nprocesses=2
nthreads=2

# Point this to your ONETEP executable 
# (normally found in your installation's ./bin directory).
onetep_exe=/home/user/onetep/bin/onetep.gfortran

# Point this to your onetep_launcher
# (normally found in your installation's ./utils directory).
onetep_launcher=/home/user/onetep/utils/onetep_launcher

# Load any necessary modules here, set up any necessary environment variables.
# e.g. for Intel you might have to
# module load oneapi
# source /local/software/intel/oneapi/2023.0/setvars.sh

# This runs ONETEP. On some systems you may need to replace 'mpirun'
# with a suitable alternative.
mpirun -np 2 $onetep_launcher -t $nthreads -e $onetep_exe $@