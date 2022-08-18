#!/bin/bash
# "=" should not surround by space in variable assignment
l=6

fpatten="/data/inspur_disk01/userdir/jyj/FAST/jingyj/G15_new/G15_${l}/data2/G15_${l}_arcdrift-M*_W-2021*-specs_T.hdf5"

fpaths="$(ls $fpatten)"
# echo "$fpaths"

outdir="/data/inspur_disk01/userdir/xuc/FAST/G15/new_sep/drift${l}/"
outdirs=${outdir}"sub"
outdirr=${outdir}"rfi"
outdirf=${outdir}"fft"
outdirm=${outdir}"mul_fft"
for dir in ${outdirs} ${outdirr} ${outdirf} ${outdirm} 
    do
        if [ ! -d $dir ];then
           mkdir -p $dir
           echo "makedir $dir"
        fi
    done

commands=$(cat <<EOF
python -m hifast.bld | -c ../conf/G15-1_bld.ini -f --outdir $outdirs
python -m hifast.bld | -c ../conf/G15-med_2_bld.ini -f
#
python -m hifast.rfi | -c ../conf/G15-2_rfi_d6.ini -f --outdir $outdirr
# 
python -m hifast.sw | -c ../conf/G15-3_fft_sw.ini -f --outdir $outdirf
python -m hifast.bld | -c ../conf/G15-1_bld.ini -f 
python -m hifast.bld | -c ../conf/G15-med_2_bld.ini -f
# #
python -m hifast.sw | -c ../conf/G15-3_fft_sw2.ini -f --outdir $outdirm
python -m hifast.flux | -c ../conf/G15-4_flux.ini -f 
python -m hifast.multi | -c ../conf/G15-5_multi.ini -f
python -m hifast.downsample | -c ../conf/G15-6_down.ini -f 
EOF
)
# echo "$commands"

hifast.sh "${fpaths}" -c "$commands" -n 2

