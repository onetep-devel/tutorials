#!/bin/bash
module load oneapi
source /local/software/intel/oneapi/2023.0/setvars.sh
ulimit -s unlimited
export OMP_NUM_THREADS=2
mpirun -np 2 /local/disk2/tom_stuff_please_do_not_delete/onetep/bin/onetep.RH8.ifort.omp.scalapack $@