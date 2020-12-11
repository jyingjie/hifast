#!/bin/bash

export OPENBLAS_NUM_THREADS=10
export MKL_NUM_THREADS=10

#for fname in  data/*M01*specs_T*.npy
for fname in  data/*M01*specs_T*.hdf5
do
if [[ $fname == *"radec"* ]]; then continue; fi
echo $fname
python -m fast_python.ky2radec $fname
done