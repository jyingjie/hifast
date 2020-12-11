#!/bin/bash

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

#for fname in data/*sub_baseline-flux*.hdf5
#for fname in data/*flux-sub_baseline*.hdf5
for fname in data/*sub_baseline*.hdf5
do
if [[ $fname == *"corr_vel"* ]]; then continue; fi
echo $fname
python -m fast_python.corr_vel $fname
done

