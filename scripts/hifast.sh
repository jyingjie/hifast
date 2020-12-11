#!/bin/bash

set -e
trap 'kill $BGPID; exit' INT


POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -n)
    nproc="$2"
    shift # past argument
    shift # past value
    ;;
    -i|--files)
    files="$2"
    shift # past argument
    shift # past value
    ;;
    -c|--command)
    command="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done

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
 if [[ ($fname != *.hdf5) && ($fname != *.fits) ]]; then echo "not hdf5 or fits file" && exit 1; fi
 if [[ ! -f $fname ]]; then echo "file ${fname} not exists" && exit 1; fi
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
sleep 1
#echo "${commands[@]}"
#
Run_fname() {
for line in "${commands[@]}"
do
  if [[ "$line" == "#"* ]]; then continue ; fi
  run_str="$(echo $line | awk -F "|" -v a="$fname" '{print($1,a,$2)}')"
  echo "run: ${run_str}"
  output="$(${run_str})"
  echo "$output"
  echo ""
  fname="$(echo "$output" | sed -En 's/Saved to |File exists //p')"
done
}

# run all files

for ((i = 0 ; i < ${#fname_list[@]} ; i++)); do
  fname="${fname_list[i]}"
  #echo $fname
  Run_fname & 
  if [[ $(( (i+1) % nproc )) -eq 0 ]]; then wait; fi
done
wait
#Run_fname &
#wc ${file_l ist[@]}
#echo $($command)

#echo $@
#echo $1