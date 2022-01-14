#!/bin/bash
# "=" should not surround by space in variable assignment
l=6

fpatten="/data/inspur_disk01/userdir/xuc/FAST/G15/new_22/drift${l}/sep/G15_${l}_arcdrift-M01_W-202108*-specs_T.hdf5"
fpaths="$(ls $fpatten)"
# echo "$fpaths"
outdir="/data/inspur_disk01/userdir/xuc/FAST/G15/new_22/drift${l}/med"
if [ ! -d $outdir ];then
   mkdir -p $outdir
fi

commands=$(cat <<EOF
python /data/inspur_disk01/userdir/xuc/FAST/test_fast/hifast_dev/hifast/sw.py | -c ../conf/G15-med_1_sw.ini \
        --outdir $outdir -f
python -m hifast.bld | -c ../conf/G15-med_2_bld.ini -f
EOF
)
# echo "$commands"

/data/inspur_disk01/userdir/xuc/FAST/G15/new_22/scripts/hifast.sh "${fpaths}" -c "$commands" -n 2

