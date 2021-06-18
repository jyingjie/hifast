#!/bin/bash

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

l=6
echo "dealing with folder drift_${l}"

: << EOF
outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/newG15/drift_${l}/sub_baseline_once"
mkdir $outdir

## first, polyfit baseline
for fname in drift_${l}/sep/G15_drift_${l}_arcdrift-M*-specs_T.hdf5
do
if [[ $fname == *"radec"* ]]; then continue; fi
echo $fname
python -m hifast.cli_baseline $fname --frange 1320 1440  \
    --nproc 10 --method arPLS --s_method_freq gaussian --s_sigma_freq 3 --lam 1e7   -f \
    --outdir $outdir -f|| exit 1
done
echo "finish once sub wave"


## second, mark RFI

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/newG15/drift_${l}/cor_vel_once"
mkdir $outdir
for fname in drift_${l}/sub_baseline_once/*-bld.hdf5
do
echo $fname
python -m hifast.cli_multi $fname --tr --tr_method smooth --tr_s_sigma_t 10 --tr_times 5 \
        --pr_times_s 1.6 --tr_n_continue 80 --ext_add 3  --cross_frac .5 \
        --flux --fc --frame HELIOCENT -f \
        --outdir  $outdir  || exit 1
done


EOF
##########################

## or second, mark RFI
outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/newG15/drift_${l}/cor_vel_once"
mkdir $outdir
for fname in drift_${l}/sub_baseline_once/*-bld.hdf5
do
echo $fname
python cli_markRFI.py $fname --outdir  $outdir \
        --time_rfi --lf_beams ['05','06','13'] \
        --sf_frange 1380 1382 --sf_times 10 --sf_thr 10 --sf_rfi_last 20 --sf_T_thr .3 \
        --lf_frange 1400 1450 --lf_times 1.5 --lf_thr 0 --lf_rfi_last 50 --lf_ext 10 \
        \
        --s_method_freq gaussian --s_sigma_freq 3 --s_method_t boxcar --s_sigma_t 7 \
        --rfi_thr 3 --rms_frange 1400 1403 --mw_frange 1419 1425 --rfi_groups 'two groups' \
        --rfi_width_lim 20 --ext_sec 20 --freq_thr .3 --freq_step 8.1 \
        --ext_edge 20 --mask_thr 15 \
        --plot -f  || exit 1 
done
## third, remove ripple

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/newG15/drift_${l}/sub_baseline_once"
#rfi_freq_step=1/16.2

for fname in drift_${l}/sub_baseline_once/*-bld.hdf5
do
echo $fname
python cli_baseline_fft.py $fname --outdir $outdir \
        --fft_method rfft --sw_freq 0.9254 --rfi_freq_step 0.0617283950617284  --amp_thr 35  --sw_n 5 \
        --rfi_method subtract --mw_frange 1420.2 1420.55 --sg_window 1.0 --sg_polyorder 7 --mw_lower 1.0e4 \
        --plot -f || exit 1 
done
echo "finish standing wave subtraction"



## fourth, flux and velocity correction

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/newG15/drift_${l}/cor_vel_once"
#mkdir $outdir
for fname in drift_${l}/sub_baseline_once/*-fft_rfi_bld.hdf5
do
echo $fname
python -m hifast.cli_multi $fname  --flux --fc --frame HELIOCENT -f \
        --outdir  $outdir --keep_polar || exit 1
done
        
echo "finish flux and velocity correction"


outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/newG15/drift_${l}/cor_vel_once/down"
mkdir $outdir
for fname in drift_$l/cor_vel_once/*fft*flux_fc.hdf5

do
echo $fname
python  downsample.py $fname \
        --outdir $outdir \
        --chanfactor 3 -f || exit 1
done
echo "finish down sample"

