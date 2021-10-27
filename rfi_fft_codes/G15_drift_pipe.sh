#!/bin/bash

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

folder_num=7
echo "dealing with folder drift ${folder_num}"

## after sep ...
: << EOF
outdir="/data/inspur_disk01/userdir/xuc/FAST/G15/new_21/drift${folder_num}/sub"
mkdir $outdir

## first, polyfit baseline
for fname in /data/inspur_disk01/userdir/jyj/FAST/jingyj/G15_new/G15_${folder_num}/data/G15_*_arcdrift-M*_W-2021*-specs_T.hdf5
do
if [[ $fname == *"radec"* ]]; then continue; fi
echo $fname
python -m hifast.cli_baseline $fname --frange 1300 1460  \
    --s_method_freq gaussian --s_sigma_freq 3 \
    --nproc 10 --method arPLS  --lam 1e9 -f \
    --outdir $outdir -f|| exit 1
done

echo "1. finish once sub wave"

#############################

outdir="/data/inspur_disk01/userdir/xuc/FAST/G15/new_21/drift${folder_num}/rfi"
mkdir $outdir

for fname in ./drift${folder_num}/sub/G15_*_arcdrift-M*_W-2021*-specs_T-bld.hdf5
do
echo $fname
python /data/inspur_disk01/userdir/xuc/FAST/test_fast/test_pipe/cli_mark_tRFI.py $fname \
        --mw_frange 1420.2 1420.45 --lf_sepname 'input_subname' \
        --lf --lf_frange 1385 1460 --lf_times 2 --lf_thr 0 --lf_rfi_last 50 --lf_ext 50 \
        --sf --sf_frange_step 30 --sf_times 3 --sf_thr 1 --sf_rfi_last 10 --sf_T_thr_times 2.5 \
        --plot --vmin_max -.05 .05 -f  --outdir $outdir || exit 1 
done

echo "2. mark time rfi"
#EOF
#############################

outdir="/data/inspur_disk01/userdir/xuc/FAST/G15/new_21/drift${folder_num}/sub_fft"
mkdir $outdir
sepname="/data/inspur_disk01/userdir/jyj/FAST/jingyj/G15_new/G15_${folder_num}/data/G15_${folder_num}_arcdrift-M"

for fname in ./drift${folder_num}/sub/G15_*_arcdrift-M*_W-2021*-specs_T-bld.hdf5
do
echo $fname
python /data/inspur_disk01/userdir/xuc/FAST/test_fast/test_pipe/cli_baseline_fft.py $fname \
        -sep $sepname  --outdir $outdir  \
        --mw_frange  1420.2 1420.45 --rms_sigma 6 --rms_frange 1390 1400 \
        --rfi_method 'near ripple' --times_lower_thr 4.5 \
        --rfi_width_lim 10 --ext_sec 5 --ext_freq 1.3 \
        --fft_method rfft --chan_wide 7 --chan_narr 3 \
        --amp_thr_mean_factor 1.1 --amp_thr_factor 1.4 --choose_method 'all' \
        --rip_base --rip_1mhz --rip_2mhz --rip_0_04mhz --fft_ylim -5 160 \
        --plot --ylim -1 .5 --vmin_max -.05 .05 -f --fill_rfi rfi --keep_polar || exit 1 
done

echo "3. finish subtracting fft ripple"
EOF

###########################

outdir="/data/inspur_disk01/userdir/xuc/FAST/G15/new_21/drift${folder_num}/sub_fft"
for fname in ./drift${folder_num}/sub_fft/G15_*_arcdrift-M*_W-2021*-specs_T-bld-fft_bldea.hdf5
do
echo $fname
python -m hifast.cli_baseline $fname \
    --s_method_freq gaussian --s_sigma_freq 3 \
    --nproc 10 --method arPLS  --lam 1e9 -f \
    --outdir $outdir -f|| exit 1

done
echo "4. finish sub again"


####################################
outdir="/data/inspur_disk01/userdir/xuc/FAST/G15/new_21/drift${folder_num}/sub_fft"

for fname in ./drift${folder_num}/sub_fft/G15_*_arcdrift-M*_W-2021*-specs_T-bld-fft_bldea-bld.hdf5
do
echo $fname
python /data/inspur_disk01/userdir/xuc/FAST/test_fast/test_pipe/cli_baseline_fft.py $fname \
        -sep 'none'  --outdir $outdir --frange 1300 1439 \
        --mw_frange  1420.2 1420.45 --rms_sigma 6 --rms_frange 1390 1400 \
        --rfi_method 'near ripple' --times_lower_thr 4.5 \
        --rfi_width_lim 10 --ext_sec 5 --ext_freq 1.3 \
        --fft_method rfft \
        --amp_thr_mean_factor 1.1 --amp_thr_factor 1.4 --choose_method 'all' \
        --rip_base --fft_ylim -5 160 \
        --plot --ylim -1 .5 --vmin_max -.05 .05 -f --fill_rfi nan --keep_polar || exit 1 
done

echo "4.5. finish subtracting fft baseline"
: << EOF
###########################


outdir="/data/inspur_disk01/userdir/xuc/FAST/G15/new_21/drift${folder_num}/cor_vel"
mkdir $outdir
for fname in ./drift${folder_num}/sub_fft/G15_*_arcdrift-M*_W-2021*-specs_T-bld-fft_bldea-bld-fft_bldea.hdf5
do
echo $fname
python -m hifast.cli_multi $fname  --flux --fc --frame HELIOCEN \
        --pr --pr_times_s 1.5 --ext_frac .3 \
        --outdir  $outdir  -f || exit 1
done       
        
echo "5. finish flux and velocity correction, polar rfi"


#######################

outdir="/data/inspur_disk01/userdir/xuc/FAST/G15/new_21/drift${folder_num}/cor_vel/down"
mkdir $outdir
for fname in ./drift${folder_num}/cor_vel/G15_*_arcdrift-M*_W-2021*-specs_T-bld-fft_bldea-bld-fft_bldea-rfi_flux_fc.hdf5
do
echo $fname
python  /data/inspur_disk01/userdir/xuc/FAST/test_fast/test_pipe/downsample.py $fname --outdir $outdir \
        --chanfactor 3 -f || exit 1
done

echo "6. finish down sample"
EOF
