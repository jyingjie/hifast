#!/bin/bash

fpart="/data/hw1/FAST/3047/M31_Drift_v3_12/20201201/M31_Drift_v3_12_arcdrift-M"

for i in {01..02}
do
python -m hifast.sep "${fpart}${i}_W_0001.fits" -d 6 -m 6 -n 1994 --step 5 --frange 1414 1431 \
                                --smooth poly --s_deg 1 --outdir ./data_v3 --save_p_cal -f
done
