#!/bin/bash

export OPENBLAS_NUM_THREADS=10
export MKL_NUM_THREADS=10
#python -m fast_python.fits_cubes --help
python -m fast_python.fits_cubes **/data/*-corr_vel*.hdf5 \
                          --outname ./test_cubes.fits \
                          --bwidth 60 -p SIN