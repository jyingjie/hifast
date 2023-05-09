#!/usr/bin/env bash

set -e
trap 'kill $BGPID; exit' INT


POSITIONAL=()
while [[ $# -gt 0 ]]
do
key=$1
# https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
case $key in
    -n)
    nproc=$2
    shift # past argument
    shift # past value
    ;;
    -i|--files)
    files=$2
    shift # past argument
    shift # past value
    ;;
    -c|--command)
    command=$2
    shift # past argument
    shift # past value
    ;;
    -s|--save_log)
    save_log=YES
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=($1) # save it in an array for later
    shift # past argument
    ;;
esac
done

echo "-s = $save_log"
#set -- "${POSITIONAL[@]}" # restore positional parameters

# set default nproc
if [[ -z $nproc ]]; then
nproc=1
fi
#echo $nproc
# load fnames

if [[ -n $files ]]; then
 fname_list=()
 if [[ ! -f $files ]]; then echo "input file lists not exists" && exit 1; fi
 while IFS= read -r line
 do
 fname_list+=("$line")
 done < "$files"
else
 fname_list=("${POSITIONAL[@]}")
fi

#checking
echo "processing:"
for fname in "${fname_list[@]}"
do
# if [[ ($fname != *.hdf5) && ($fname != *.fits) ]]; then echo "not hdf5 or fits file" && exit 1; fi
 if [[ ! -f $fname ]] && [[ ! -d $fname ]]; then echo "file ${fname} not exists" && exit 1; fi
 echo "$fname"
done
sleep 1
#commands
echo "Commands:"
commands=()
if [[ "$command" == *".par" ]]; then
  echo "load commands from ${command}"
  while IFS= read line || [ -n "$line" ]
    do
    if [[ "$line" == "#"* ]]; then continue ; fi
    echo "$line"
    commands+=("$line")
    done < "$command"
elif [[ "$command" == *"|"* ]]; then
  commands+=("$command")
fi
coms="$(printf "%s\n" "${commands[@]}")"

sleep 1
#echo "${commands[@]}"
#
Run_fname() {
local fname="$1"
echo "$coms" | while read line; do
  if [[ "$line" == "#"* ]]; then continue ; fi
  run_str="$(echo "$line" | awk -F "|" -v a="$fname" '{print($1,a,$2)}')"
  echo "run: ${run_str}"
  output="$(${run_str})" || exit 1
  echo "$output"
  echo ""
  fname="$(echo "$output" | sed -En 's/Saved to |File exists //p')"
done
}

# run all files
export -f Run_fname
export coms
if [[ $save_log == "YES" ]]; then
  lognameadd=`date +"%Y%m%d-%H%M%S"`
  export lognameadd
  mkdir -p log
  echo "will save output to directory log"
  echo "${fname_list[@]}" | xargs -n 1 -P $nproc bash -c \
    'echo "$coms"; Run_fname "$1" > >( tee -a "log/$(basename $1).$lognameadd.out") 2> >(tee >(grep -v '███' >> "log/$(basename $1).$lognameadd.err") >&2)' _
else
  echo "${fname_list[@]}" | xargs -n 1 -P $nproc bash -c 'echo "$coms"; Run_fname "$1"' _
fi
#Run_fname &
#wc ${file_l ist[@]}
#echo $($command)

#echo $@
#echo $1
