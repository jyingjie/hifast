#!/bin/bash

export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1


echo "dealing with folder drift 6"

## after sep ...
#: << EOF
outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/new_21/drift6/sub"
mkdir $outdir

## first, polyfit baseline
for fname in /data/inspur_disk01/userdir/jyj/FAST/jingyj/G15_new/G15_6/data/G15_*_arcdrift-M*_W-2021*-specs_T.hdf5
do
if [[ $fname == *"radec"* ]]; then continue; fi
echo $fname
python -m hifast.cli_baseline $fname --frange 1300 1460  \
    --nproc 10 --method arPLS  --lam 1e7 -f \
    --outdir $outdir -f|| exit 1
done

echo "1. finish once sub wave"


#############################

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/new_21/drift6/rfi"
mkdir $outdir
sepname='/data/inspur_disk01/userdir/jyj/FAST/jingyj/G15_new/G15_6/data/G15_6_arcdrift-M'
for fname in ./drift6/sub/G15_*_arcdrift-M*_W-2021*-specs_T-bld.hdf5
do
echo $fname
python /data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/test_pipe/cli_mark_tRFI.py $fname \
        --mw_frange 1420.2 1420.45 --lf_sepname $sepname \
        --lf --lf_frange 1421 1450 --lf_times 1.05 --lf_thr 0 --lf_rfi_last 500 --lf_ext 1 \
        --plot --vmin_max -.05 .05 -f  --outdir $outdir || exit 1 
done

echo "2. mark time rfi"

#############################

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/new_21/drift6/sub_fft"
mkdir $outdir
sepname='/data/inspur_disk01/userdir/jyj/FAST/jingyj/G15_new/G15_6/data/G15_6_arcdrift-M'

for fname in ./drift6/sub/G15_*_arcdrift-M*_W-2021*-specs_T-bld.hdf5
do
echo $fname
python /data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/test_pipe/cli_baseline_fft.py $fname \
        -sep $sepname  --outdir $outdir  \
        --mw_frange  1420.2 1420.45 --rms_sigma 6 --rms_frange 1390 1400 \
        --rfi_method 'near ripple' --times_lower_thr 4 \
        --rfi_width_lim 15 --ext_sec 20 --ext_freq 1.3 \
        --fft_method rfft --chan_wide 7 --chan_narr 3 \
        --amp_thr_mean_factor 1.05 --amp_thr_factor 1.4 --choose_method 'all' \
        --rip_base --rip_1mhz --rip_2mhz --rip_0_04mhz --fft_ylim -5 160 \
        --plot --ylim -1 .5 --vmin_max -.05 .05 -f --fill_rfi nan --keep_polar || exit 1 
done

echo "3. finish subtracting fft ripple"

###########################

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/new_21/drift6/sub_fft"
for fname in ./drift6/sub_fft/G15_*_arcdrift-M*_W-2021*-specs_T-bld-fft_bldea.hdf5
do
echo $fname
python -m hifast.cli_baseline $fname \
    --s_method_freq gaussian --s_sigma_freq 3 \
    --nproc 10 --method arPLS  --lam 1e7 -f \
    --outdir $outdir -f|| exit 1

done
echo "4. finish sub again"


####################################
outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/new_21/drift6/sub_fft"

for fname in ./drift6/sub_fft/G15_*_arcdrift-M*_W-2021*-specs_T-bld-fft_bldea-bld.hdf5
do
echo $fname
python /data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/test_pipe/cli_baseline_fft.py $fname \
        -sep 'none'  --outdir $outdir  \
        --mw_frange  1420.2 1420.45 --rms_sigma 6 --rms_frange 1390 1400 \
        --rfi_method 'near ripple' --times_lower_thr 4 \
        --rfi_width_lim 15 --ext_sec 20 --ext_freq 1.3 \
        --fft_method rfft \
        --amp_thr_mean_factor 1.05 --amp_thr_factor 1.4 --choose_method 'all' \
        --rip_base --fft_ylim -5 160 \
        --plot --ylim -1 .5 --vmin_max -.05 .05 -f --fill_rfi nan --keep_polar || exit 1 
done

echo "4.5. finish subtracting fft baseline"
# EOF
###########################


outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/new_21/drift6/cor_vel"
mkdir $outdir
for fname in ./drift6/sub_fft/G15_*_arcdrift-M*_W-2021*-specs_T-bld-fft_bldea-bld-fft_bldea.hdf5
do
echo $fname
python -m hifast.cli_multi $fname  --flux --fc --frame HELIOCEN \
         \
        --pr --pr_times_s 2 --ext_add 1 \
        --outdir  $outdir  -f  || exit 1
done
        
echo "5. finish flux and velocity correction, polar rfi"

#######################

outdir="/data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/new_21/drift6/cor_vel/down"
mkdir $outdir
for fname in ./drift6/cor_vel/G15_*_arcdrift-M*_W-2021*-specs_T-bld-fft_bldea-bld-fft_bldea-rfi_flux_fc.hdf5
do
echo $fname
python  /data/inspur_disk01/userdir/ucas_students/xuc/FAST/G15/test_pipe/downsample.py $fname --outdir $outdir \
        --chanfactor 3 -f || exit 1
done

echo "6. finish down sample"
