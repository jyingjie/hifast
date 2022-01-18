#!/bin/bash

export OPENBLAS_NUM_THREADS=10
export MKL_NUM_THREADS=10

outdir="/data/inspur_disk01/userdir/xuc/FAST/G15/new_22/cubes/"
if [ ! -d $outdir ];then
   mkdir -p $outdir
fi

python -m hifast.cube /data/inspur_disk01/userdir/xuc/FAST/G15/new_22/drift6/mul_fft/down/G15*flux*.hdf5 \
                          --outname ${outdir}'G15_6_arcdrift_HEL_fft_cube.fits' \
                          --bwidth 60 \
                          -p SIN \
                          --method 'reweight' \
                          --r_cut 60 || exit 1


