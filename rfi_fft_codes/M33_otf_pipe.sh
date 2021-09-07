#!/bin/bash

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1


echo "dealing with folder OTF"

## after sep ...
#: << EOF
outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/M33/new_21/OTF_h/sub"
mkdir $outdir

## first, polyfit baseline
for fname in /data/inspur_disk01/userdir/jyj/FAST/jingyj/M33_new/data/M33_OTF_1_MultiBeamOTF-M*_W-20210731-specs_T.hdf5
do
if [[ $fname == *"radec"* ]]; then continue; fi
echo $fname
python -m hifast.cli_baseline $fname --frange 1300 1460  \
    --nproc 10 --method arPLS  --lam 1e7 -f \
    --outdir $outdir -f|| exit 1
done

echo "1. finish once sub wave"


#############################

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/M33/new_21/OTF_h/sub_fft"
mkdir $outdir
sepname='/data/inspur_disk01/userdir/jyj/FAST/jingyj/M33_new/data/M33_OTF_1_MultiBeamOTF-M'

for fname in ./OTF_h/sub/M33_OTF_1_MultiBeamOTF-M*_W-20210731-specs_T-bld.hdf5
do
echo $fname
python /data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/test_pipe/cli_baseline_fft.py $fname \
        -sep $sepname   -rfi 'none' --outdir $outdir  \
        --mw_frange 1420.2 1420.55 --rms_sigma 6 --rms_frange 1390 1400 \
        --rfi_method 'near ripple' --times_lower_thr 3.5 \
        --rfi_width_lim 15 --ext_sec 20 --ext_freq 1.3 \
        --fft_method rfft --chan_wide 5 --chan_narr 3 \
        --amp_thr_mean_factor 1.05 --amp_thr_factor 1.4 --choose_method 'all' \
        --rip_base --rip_1mhz --rip_2mhz --rip_0_04mhz --fft_ylim -5 140 \
        --plot --ylim -1 .5 --vmin_max -.05 .05 -f --fill_rfi nan --keep_polar || exit 1 

done

echo "2. finish subtracting fft ripple"

###########################

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/M33/new_21/OTF_h/sub_fft"
for fname in ./OTF_h/sub_fft/M33_OTF_1_MultiBeamOTF-M*_W-20210731-specs_T-bld-fft_bldea.hdf5
do
echo $fname
python -m hifast.cli_baseline $fname \
    --s_method_freq gaussian --s_sigma_freq 3 \
    --nproc 10 --method arPLS  --lam 1e7 -f \
    --outdir $outdir -f|| exit 1

done
echo "3. finish sub again"
#EOF

############################

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/M33/new_21/OTF_h/cor_vel"
mkdir $outdir
for fname in ./OTF_h/sub_fft/M33_OTF_1_MultiBeamOTF-M*_W-20210731-specs_T-bld-fft_bldea-bld.hdf5
do
echo $fname
python -m hifast.cli_multi $fname  --fc --frame LSRK \
        --pr --pr_times_s 2 --ext_add 1 \
        --outdir  $outdir  -f  || exit 1
# --flux -c '/data/inspur_disk01/userdir/jyj/FAST/jingyj/M33_new/3C48_OTF_MultiBeamOTF_20210731_K2Jy_Gain.hdf5' \
done
        
echo "4. finish flux and velocity correction, polar rfi"

#######################
: << EOF
outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/M33/new_21/OTF_h/cor_vel/down"
mkdir $outdir
for fname in ./OTF_h/cor_vel/M33_OTF_1_MultiBeamOTF-M*_W-20210731-specs_T-bld-fft_bldea-bld-rfi_flux_fc.hdf5
do
echo $fname
python  ../../G15/newG15/downsample.py $fname --outdir $outdir \
        --chanfactor 3 -f || exit 1
done

echo "5. finish down sample"
#EOF
