#!/bin/bash
# "=" should not surround by space in variable assignment
l=6

# fpatten="/data/inspur_disk01/userdir/xuc/FAST/G15/new_22/drift${l}/sep/G15_${l}_arcdrift-M*_W-202108*-specs_T.hdf5"
 fpatten="/data/inspur_disk01/userdir/xuc/FAST/G15/new_22/drift${l}/sub/G15_${l}_arcdrift-M*_W-202108*-specs_T-bld-bld.hdf5"
fpaths="$(ls $fpatten)"
# echo "$fpaths"
routdir="/data/inspur_disk01/userdir/xuc/FAST/G15/new_22/drift${l}/rfi"
if [ ! -d $routdir ];then
   mkdir -p $routdir
fi

commands=$(cat <<EOF
# python -m hifast.bld | -c ../conf/G15-1_bld.ini --outdir $outdir -f
# python -m hifast.bld | -c ../conf/G15-med_2_bld.ini -f
python -m hifast.rfi | -c ../conf/G15-2_rfi.ini --outdir $routdir -f
EOF
)
# echo "$commands"

/data/inspur_disk01/userdir/xuc/FAST/G15/new_22/scripts/hifast.sh "${fpaths}" -c "$commands" -n 2

