#!/usr/bin/env bash

# Strict mode
set -euo pipefail
IFS=$'\n\t'
trap 'kill $(jobs -p) 2>/dev/null || true' EXIT INT TERM

# Global variables
declare -r VERSION="0.0.0"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAX_PROCS=$(nproc)
DEFAULT_NPROC=1
VERBOSE=0
export HIFAST_MAKE_OUTPUT_CLEAR=${HIFAST_MAKE_OUTPUT_CLEAR:-1}

# Error handling function
error() {
    echo "ERROR: $1" >&2
    exit "${2:-1}"
}

# Warning function
warn() {
    echo "WARNING: $1" >&2
}

# Debug function
debug() {
    if ((VERBOSE)); then
        echo "DEBUG: $1" >&2
    fi
}

# Help message function (previous implementation)
show_help() {
    cat << EOF
$(basename "$0") - Execute commands in parallel on multiple files.

Usage: 
    $(basename "$0") [options] [files...]
    $(basename "$0") -i <file_list_file> [options]
    $(basename "$0") -h | --help

Description:
    Execute commands in parallel on multiple files, primarily designed for 
    processing FAST data files. Supports both direct command execution 
    and parameter file-based processing.

Options:
    -h, --help             Show this help message
    -n <number>            Number of parallel processes (default: 1)
    -i, --files <file>     File containing list of input files
    -c, --command <cmd>    Commands in a variable or parameter file (.par) to execute
    -s, --save_log         Save output to log files (in ./log directory)

Input Methods for files:
    1. Direct input:       Provide files as arguments
    2. File list:          Use -i option with a file containing paths (one per line)

Command Formats for command (-c):
    1. Parameter file (.par):
       - Each line is executed as a separate command
       - Lines starting with # are treated as comments
       - The input filename is automatically inserted between (|) symbols
       - Example: my_commands.par:
         # Process data
         python process.py | --params 
         python analyze.py | --other-params
    
    2. Commands in a variable:
       - Example:
        commands=\$(cat << 'EOF'
        python process.py | --params 
        python analyze.py | --other-params
        EOF
        )
       

Examples:
    # Basic usage with parameter file
    $(basename "$0") -n 3 -c params.par file1.hdf5 file2.hdf5

    # Process HDF5 files with direct command and logging
    $(basename "$0") -c "python -m hifast.cbr | --outdir ./%[date]s" -s *.hdf5

    # Process files from a list with 4 parallel processes
    $(basename "$0") -i file_list.txt -c command.par -n 4

    # Real-world example for FAST data processing
    files=\$(ls /path/to/FAST/data/*.hdf5)
    commands=\$(cat << 'EOF'
    python -m hifast.cbr | --outdir ./%[project]s/%[date]s \\
        --nproc 9 --plot --obsmode MultiBeamCalibration \\
        --cbrname SOURCE --frange 1290 1450
    EOF
    )
    $(basename "$0") "\$files" -c "\$commands" -n 1 -s

Log Files:
    When -s is used, logs are saved to:
        ./log/<filename>.<timestamp>.out  - Standard output
        ./log/<filename>.<timestamp>.err  - Error messages
    Timestamp format: YYYYMMDD-HHMMSS

Common Issues:
    1. For large files, adjust -n based on available system memory
    2. Use quotes around commands containing special characters

Environment:
    The script respects the following environment variables:
    - TMPDIR: Temporary directory for intermediate files
    - PATH: For finding executables in commands

See Also:
    hifast documentation: https://hifast.readthedocs.io/


EOF
    exit 0
}

# Function to validate number input
validate_number() {
    local num=$1
    local name=$2
    if ! [[ "$num" =~ ^[0-9]+$ ]]; then
        error "$name must be a positive integer" 2
    fi
}

# Function to validate file existence
validate_file() {
    local file=$1
    if [[ ! -f "$file" ]]; then
        error "File not found: $file" 2
    fi
    if [[ ! -r "$file" ]]; then
        error "File not readable: $file" 2
    fi
}

Run_fname() {
    local fname="$1"
    local start_time
    start_time=$(date +%s)
    local status=0
    local current_fname="$fname"
    
    debug "Processing file: $fname"
    
    # Process each command line
    while IFS= read line || [ -n "$line" ]; do
        # Skip comments
        [[ "$line" =~ ^[[:space:]]*# ]] && continue
        
        # Replace the pipe symbol with the filename
        if [[ "$line" == *"|"* ]]; then
            # Split command into parts and construct final command
            local cmd_part1 cmd_part2
            cmd_part1=$(echo "$line" | cut -d'|' -f1)
            cmd_part2=$(echo "$line" | cut -d'|' -f2-)
            
            # Construct command without eval
            set -- $cmd_part1 "$current_fname" $cmd_part2
            
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running: $*"
            
            # Try using fd3 with error handling
            if exec 3>&1; then
                output=$("$@" | tee /dev/fd/3)
                exit_code=${PIPESTATUS[0]}
                echo "" >&3
                exec 3>&-  # Close fd3
            else
                output=$("$@")
                exit_code=$?
                echo "$output"
                echo ""
            fi
            
            if [ $exit_code -ne 0 ]; then
                printf "Command failed:\n%s\n" "$*" >&2
                exit $exit_code
            fi
            
            # Update filename for next command if needed
            if [[ "$output" =~ Saved[[:space:]]to[[:space:]]([^[:space:]]+) ]]; then
                current_fname="${BASH_REMATCH[1]}"
            elif [[ "$output" =~ File[[:space:]]exists[[:space:]]([^[:space:]]+) ]]; then
                current_fname="${BASH_REMATCH[1]}"
            else
                printf "No file found for output. Command failed:\n%s\n" "$*" >&2
                exit 1
            fi
        fi
    done < <(if [[ "$command" == *.par ]]; then cat "$command"; else echo "$command"; fi)
    
    local end_time
    end_time=$(date +%s)
    debug "Processing time for $fname: $((end_time - start_time)) seconds"
    return $status
}


# Parse command line arguments
POSITIONAL=()
nproc=$DEFAULT_NPROC
save_log=""
files=""
command=""

# Show help if no parameters
if [ $# -eq 0 ]; then
    show_help
fi

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            ;;
        -v|--verbose)
            VERBOSE=1
            shift
            ;;
        -n)
            nproc=$2
            validate_number "$nproc" "Number of processes"
            shift 2
            ;;
        -i|--files)
            files=$2
            validate_file "$files"
            shift 2
            ;;
        -c|--command)
            command=$2
            [[ -n $command ]] || error "Command cannot be empty"
            [[ "$command" == *.par ]] && validate_file "$command"
            shift 2
            ;;
        -s|--save_log)
            save_log=YES
            mkdir -p log || error "Cannot create log directory"
            shift
            ;;
        *)
            POSITIONAL+=($1)
            shift
            ;;
    esac
done

# Validate inputs
if [[ ${#POSITIONAL[@]} -eq 0 ]] && [[ -z $files ]]; then
    error "No input files specified"
fi

# Load file list
if [[ -n $files ]]; then
    mapfile -t fname_list < "$files"
else
    fname_list=("${POSITIONAL[@]}")
fi

# Validate all files before processing
for fname in "${fname_list[@]}"; do
    validate_file "$fname"
done

# Export necessary functions and variables for xargs subshell
export -f Run_fname debug
export command VERBOSE

# Process files
total_files=${#fname_list[@]}
echo "###############################Files############################################"
echo "Processing $total_files files with $nproc parallel processes:"
for fname in "${fname_list[@]}"
do
# if [[ ($fname != *.hdf5) && ($fname != *.fits) ]]; then echo "not hdf5 or fits file" && exit 1; fi
 if [[ ! -f $fname ]] && [[ ! -d $fname ]]; then echo "file ${fname} not exists" && exit 1; fi
 echo "$fname"
done
echo "################################################################################"

if [[ $save_log == "YES" ]]; then
    lognameadd=$(date +"%Y%m%d-%H%M%S")
    echo "Output will be saved to directory ./log"
    echo ""
    # Process files in parallel with logging
    printf '%s\n' "${fname_list[@]}" | \
    xargs -n 1 -P "$nproc" -I {} bash -c \
        'fname="{}"; \
         exec 1> >(tee "log/$(basename "$fname").'$lognameadd'.out") \
         2> >(tee >(grep -v '███' >> "log/$(basename "$fname").'$lognameadd'.err") >&2); \
         Run_fname "$fname"'
else
    # Process files in parallel without logging
    printf '%s\n' "${fname_list[@]}" | \
    xargs -n 1 -P "$nproc" -I {} bash -c 'Run_fname "{}"'
fi

echo "Processing completed"
