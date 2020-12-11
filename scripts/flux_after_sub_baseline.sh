#!/bin/bash
for fname in ./data/*sub_baseline*.hdf5
do
  if [[ $fname == *"corr_vel"* ]]; then continue; fi
  if [[ $fname == *"flux"* ]]; then continue; fi
  echo $fname
  python -m fast_python.flux $fname
done