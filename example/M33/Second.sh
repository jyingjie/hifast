#!/bin/bash
# "=" should not surround by space in variable assignment
fpatten="./data/M33_OTF_1_*M0[1-2]*_W*-specs_T.hdf5"
fpaths="$(ls $fpatten)"
echo "$fpaths"

commands=$(cat <<EOF
python -m hifast.flux | -c ./conf/S1-flux.ini --outdir ./data_S 
python -m hifast.bld | -c ./conf/S2-bld.ini -f
python -m hifast.rfi | -c ./conf/S3-rfi.ini
python -m hifast.sw | -c ./conf/S4-sw.ini
python -m hifast.multi | -c ./conf/S5-multi.ini 
EOF
)
echo "$commands"

hifast.sh "${fpaths}" -c "$commands" -n 2
