#!/bin/bash

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

l=6
echo "dealing with folder drift_${l}"

## after sep ...
: << EOF
outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/newG15/drift_${l}/sub_baseline_fft"
mkdir $outdir

## first, polyfit baseline
for fname in drift_${l}/sep/G15_drift_${l}_arcdrift-M*-specs_T.hdf5
do
if [[ $fname == *"radec"* ]]; then continue; fi
echo $fname
python -m hifast.cli_baseline $fname --frange 1300 1440  \
    --nproc 10 --method arPLS --s_method_freq gaussian --s_sigma_freq 3 --lam 1e7   -f \
    --outdir $outdir -f|| exit 1
done

echo "1. finish once sub wave"
#EOF

#: << EOF
##########################

## second, mark RFI for fft
outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/newG15/drift_${l}/rfi"
# mkdir $outdir
for fname in drift_${l}/sub_baseline_fft/*T-bld.hdf5
do
echo $fname
python cli_markRFI.py $fname --outdir  $outdir \
        --sf --lf --lf_beams ['05','06','13'] \
        --sf_frange 1380 1382 --sf_times 10 --sf_thr 10 --sf_rfi_last 20 --sf_T_thr_times 3 \
        --lf_frange 1400 1450 --lf_times 1.5 --lf_thr 0 --lf_rfi_last 50 --lf_ext 10 \
        --period_rfi \
        --s_method_freq gaussian --s_sigma_freq 3 --s_method_t boxcar --s_sigma_t 7 \
        --rfi_thr 14 --rms_frange 1390 1400 --mw_frange 1419 1425 --rfi_groups 'two groups' \
        --rfi_width_lim 20 --ext_sec 20 --freq_thr .7 --freq_step 8.1 \
        --mask_RFI_method 'fixed freq' --ext_edge 0 --mask_thr 16 --freq_from_theory .5 \
        --plot --ylim -1 5 --save_sf -f --keep_polar --time_coherent_per .99 || exit 1 
done

echo "2. finish mark RFI"
#EOF

#: << EOF
#####################
## third, remove ripple

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/newG15/drift_${l}/sub_baseline_fft"
#rfi_freq_step=1/16.2

for fname in drift_${l}/sub_baseline_fft/*T-bld.hdf5
do
echo $fname
python cli_baseline_fft.py $fname --outdir $outdir --frange 1300 1440  \
        --mw_frange 1420.2 1420.55 --rms_sigma 6 --rms_frange 1390 1400 \
        --rfi_method 'subtract rfi' --sg_window 1.0 --sg_polyorder 7 \
        --times_lower 1.0e3 --times_lower_thr 6   \
        --fft_method rfft --sw_freq 0.9254   --amp_thr 40  --sw_n 5 \
        --plot --ylim -3 2 --one_spec -f --fill_rfi nan --keep_polar || exit 1 
done
echo "3. finish standing wave subtraction"
#EOF
#: << EOF
######################
## fourth, sub again
outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/newG15/drift_${l}/sub_baseline_fft"

for fname in drift_${l}/sub_baseline_fft/*T-bld-fft_bldr.hdf5
do
echo $fname
python -m hifast.cli_baseline $fname  \
    --nproc 10 --method arPLS --s_method_freq gaussian --s_sigma_freq 3 --lam 1e7   -f \
    --outdir $outdir -f|| exit 1
done
echo "4. finish sub again"
#EOF
#: << EOF
########################
## fifth, markRFI fully
outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/newG15/drift_${l}/rfi_full"
# mkdir $outdir
for fname in drift_${l}/sub_baseline_fft/*T-bld-fft_bldr-bld.hdf5
do
echo $fname
python cli_markRFI.py $fname --outdir  $outdir \
        --sf \
        --sf_frange 1380 1382 --sf_times 10 --sf_thr 10 --sf_rfi_last 20 --sf_T_thr_times 3 \
        --period_rfi \
        --s_method_freq gaussian --s_sigma_freq 3 --s_method_t boxcar --s_sigma_t 7 \
        --rfi_thr 16 --rms_frange 1390 1400 --mw_frange 1419 1425 --rfi_groups 'two groups' \
        --rfi_width_lim 20 --ext_sec 20 --freq_thr .7 --freq_step 8.1 \
        --mask_RFI_method 'fixed freq' --ext_edge 1 --mask_thr 5 --freq_from_theory .5 \
        --plot --ylim -1 5 -f --keep_polar --time_coherent_per .99 || exit 1 
done

echo "5. finish mark RFI"
#EOF
#: << EOF
########################
## sixth, flux and velocity correction, polar rfi

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/newG15/drift_${l}/cor_vel"
#mkdir $outdir
for fname in drift_${l}/rfi_full/*T-bld-fft_bldr-bld-tr_pdr.hdf5
do
echo $fname
python -m hifast.cli_multi $fname  --flux --fc --frame HELIOCENT\
        --pr --pr_times_s 2 --ext_add 1 \
        --outdir  $outdir  -f  || exit 1
done
        
echo "6. finish flux and velocity correction, polar rfi"
EOF
#: << EOF
#######################
## seventh, down sample
outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/newG15/drift_${l}/cor_vel/down"
mkdir $outdir
for fname in drift_$l/cor_vel/*T-bld-fft_bldr-bld-tr_pdr-rfi_flux_fc.hdf5
do
echo $fname
python  downsample.py $fname --outdir $outdir \
        --chanfactor 3 -f || exit 1
done

echo "7. finish down sample"
#EOF

## then, create cube
