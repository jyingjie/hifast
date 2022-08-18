#!/bin/bash
# "=" should not surround by space in variable assignment
l=6

fpatten="/data/inspur_disk01/userdir/xuc/FAST/G15/new_sep/drift${l}/rfi/G15_${l}_arcdrift-M*_W-202108*-specs_T-bld-bld-rfi.hdf5"
fpaths="$(ls $fpatten)"
# echo "$fpaths"
outdir="/data/inspur_disk01/userdir/xuc/FAST/G15/new_sep/drift${l}/sin"
if [ ! -d $outdir ];then
   mkdir -p $outdir
fi

commands=$(cat <<EOF
python -m hifast.sw | -c ../conf/G15-sin_1_sw.ini -f --outdir $outdir 
# # arPLS
# python -m hifast.bld | -c ../conf/G15-1_bld.ini -f 
# python -m hifast.bld | -c ../conf/G15-med_2_bld.ini -f 
# # MinMed
python -m hifast.bld | -c ../conf/G15-med_2_bld.ini -f --lam 1e11
python -m hifast.bld | -c ../conf/G15-1_bld2.ini -f --method MinMed
python -m hifast.bld | -c ../conf/G15-med_2_bld.ini -f 
#
python -m hifast.sw | -c ../conf/G15-3_fft_sw2.ini -f --rms_frange 1361 1366
python -m hifast.flux | -c ../conf/G15-4_flux.ini -f 
python -m hifast.multi | -c ../conf/G15-5_multi.ini -f
python -m hifast.downsample | -c ../conf/G15-6_down.ini -f
EOF
)
# echo "$commands"

hifast.sh "${fpaths}" -c "$commands" -n 2


