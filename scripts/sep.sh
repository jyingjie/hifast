#!/bin/bash

export OPENBLAS_NUM_THREADS=10
export MKL_NUM_THREADS=10

fpart="/data/inspur_disk06/fast_data/3047/GAMA_G15/20191215/GAMA_G15_2_arcdrift-M"

for i in {01..19}
do
python  -m fast_python.sep_spectra "${fpart}${i}_F_0001.fits" -m 10 -n 10 --step 5 --freql 1417 --freqh 1429 --smooth poly --s_deg 1 --outdir ./data
## or
#python -m fast_python.sep_spectra "${fpart}${i}_F_0001.fits" -m 10 -n 10 --step 5 --freql 1329 --freqh 1429 --smooth gaussian --s_sigma 5 --outdir ./data
## or
#python -m fast_python.sep_spectra "${fpart}${i}_F_0001.fits" -m 10 -n 10 --step 5 --freql 1329 --freqh 1429 --smooth gaussian --s_sigma 5 --dfactor 16 --med_filter --med_size 5 --outdir ./data
done