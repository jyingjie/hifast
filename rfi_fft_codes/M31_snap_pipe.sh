#!/bin/bash

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1


echo "dealing with folder snapshot"

## after sep ...
#: << EOF
outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/M31/snap_21/data6/sub"
mkdir $outdir

## first, polyfit baseline
for fname in /data/inspur_disk01/userdir/jyj/FAST/jingyj/M31_snapshot/data/M31_SnapShot_6_snapshot-M*_W-20210730-specs_T.hdf5
do
if [[ $fname == *"radec"* ]]; then continue; fi
echo $fname
python -m hifast.cli_baseline $fname --frange 1300 1460  \
    --nproc 10 --method arPLS  --lam 1e7 -f \
    --outdir $outdir -f|| exit 1
done

echo "1. finish once sub wave"
#EOF
########################

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/M31/snap_21/data6/rfi"
mkdir $outdir

## first, polyfit baseline
for fname in ./data6/sub/M31_SnapShot_6_snapshot-M*_W-20210730-specs_T-bld.hdf5
do
echo $fname
python /data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/test_pipe/cli_mark_tRFI.py $fname \
        --outdir $outdir \
        --mw_frange 1420.3 1423.3 --sf_frange_step 20 \
        --sf --sf_times 3 --sf_thr 0 --sf_rfi_last 10 --sf_T_thr_times 2.5 \
        --lf --lf_frange 1400 1450 --lf_times 1.5 --lf_thr 0 --lf_rfi_last 50 --lf_ext 10 \
        --plot --ylim -.5 .5 --vmin_max -.05 .05 -f  || exit 1 
done

echo "2. finish mark RFI"

#: << EOF
#############################

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/M31/snap_21/data6/sub_fft"
mkdir $outdir
sepname='/data/inspur_disk01/userdir/jyj/FAST/jingyj/M31_snapshot/data/M31_SnapShot_6_snapshot-M'

for fname in ./data6/sub/M31_SnapShot_6_snapshot-M*_W-20210730-specs_T-bld.hdf5
do
echo $fname
python /data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/test_pipe/cli_baseline_fft.py $fname \
        --outdir $outdir -sep $sepname  \
        --mw_frange 1420.3 1423.3 --rms_sigma 6 --rms_frange 1390 1400 \
        --rfi_method 'near ripple' --times_lower_thr 4 \
        --rfi_width_lim 15 --ext_sec 20 --ext_freq 1.3 \
        --fft_method rfft --chan_wide 5 --chan_narr 3 \
        --amp_thr_mean_factor 1.05 --amp_thr_factor 1.4 --choose_method 'all' \
        --rip_base --rip_1mhz --rip_0_04mhz --fft_ylim -5 160 \
        --plot --ylim -1 .5 --vmin_max -.05 .05 -f --fill_rfi nan --keep_polar || exit 1 

done

echo "3. finish subtracting fft ripple"
#EOF
###########################

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/M31/snap_21/data6/sub_fft"
for fname in ./data6/sub_fft/M31_SnapShot_6_snapshot-M*_W-20210730-specs_T-bld-fft_bldea.hdf5
do
echo $fname
python -m hifast.cli_baseline $fname \
    --s_method_freq gaussian --s_sigma_freq 3 \
    --nproc 10 --method arPLS  --lam 1e7 -f \
    --outdir $outdir -f|| exit 1

done
echo "4. finish sub again"

# emmmm fft again ? 


############################

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/M31/snap_21/data6/cor_vel"
mkdir $outdir
for fname in ./data6/sub_fft/M31_SnapShot_6_snapshot-M*_W-20210730-specs_T-bld-fft_bldea-bld.hdf5
do
echo $fname
python -m hifast.cli_multi $fname  --flux --fc --frame LSRK \
        \
        --pr --pr_times_s 2 --ext_add 1 \
        --outdir  $outdir  -f  || exit 1
done
        
echo "5. finish flux and velocity correction, polar rfi"
#: << EOF
#######################

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/M31/snap_21/data6/cor_vel/down"
mkdir $outdir
for fname in ./data6/cor_vel/M31_OTF_1_MultiBeamOTF-M*_W-20210731-specs_T-bld-fft_bldea-bld-rfi_flux_fc.hdf5
do
echo $fname
python  ../../G15/newG15/downsample.py $fname --outdir $outdir \
        --chanfactor 3 -f || exit 1
done

echo "6. finish down sample"
#EOF
