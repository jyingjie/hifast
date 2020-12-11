#!/bin/bash

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

for fname in data/*specs_T*.hdf5
do
if [[ $fname == *"radec"* ]]; then continue; fi
echo $fname
mpiexec -n 30 python -m mpi4py -m fast_python.sub_baseline_mpi $fname
done