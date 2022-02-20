#!/bin/bash
# "=" should not surround by space in variable assignment
l=6

fpatten="/data/inspur_disk01/userdir/xuc/FAST/G15/new_22/drift${l}/sep/G15_${l}_arcdrift-M01_W-202108*-specs_T.hdf5"
fpaths="$(ls $fpatten)"
# echo "$fpaths"

outdir="/data/inspur_disk01/userdir/xuc/FAST/G15/new_22/drift${l}/"
outdirs=${outdir}"sub"
outdirr=${outdir}"rfi"
outdirf=${outdir}"fft"
outdirm=${outdir}"mul_fft"
outdird=${outdir}"mul_fft/down"
for dir in outdirs outdirr outdirf outdirm outdird
    do
        if [ ! -d $dir ];then
           mkdir -p $dir
        fi
    done

commands=$(cat <<EOF
python -m hifast.bld | -c ../conf/G15-1_bld.ini -f --outdir $outdirs
python -m hifast.bld | -c ../conf/G15-med_2_bld.ini -f --outdir $outdirs
#
python -m hifast.rfi | -c ../conf/G15-2_rfi_d6.ini -f --outdir $outdirr
# #
python -m hifast.sw | -c ../conf/G15-3_fft_sw.ini -f --outdir $outdirf
python -m hifast.bld | -c ../conf/G15-1_bld.ini -f --outdir $outdirf
python -m hifast.bld | -c ../conf/G15-med_2_bld.ini -f --outdir $outdirf
# #
python -m hifast.sw | -c ../conf/G15-3_fft_sw2.ini -f --outdir $outdirm
python -m hifast.multi | -c ../conf/G15-4_multi.ini -f --outdir $outdirm
python -m hifast.flux | -c ../conf/G15-5_flux.ini -f --outdir $outdirm
python -m hifast.downsample | -c ../conf/G15-6_down.ini -f --outdir $outdird
EOF
)
# echo "$commands"

hifast.sh "${fpaths}" -c "$commands" -n 2

