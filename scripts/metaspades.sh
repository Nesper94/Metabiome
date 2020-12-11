#!/bin/bash
# Metagenomic assembly using metaSPAdes from SPAdes-3.12.0
# Written by: Phagomica_club

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: metabiome metaspades [Options] -i <input directory> -o <output directory> [-opts MetaSPADES_OPTIONS]"
    echo
    echo "Required:"
    echo "  -i in_dir       Input directory containing FASTQ files."
    echo "  -o out_dir      Directory in which results will be saved. This directory"
    echo "                  will be created if it doesn't exist."
    echo
    echo "Options:"
    echo "  -t NUM          Number of threads to use (default: 4)."
    echo "  -opts OPTIONS   MetaSPADES's options."
    echo "  -h, --help      Show this help"
}

# Exit if command is called with no arguments
validate_arguments "$#"

while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift ;;
        -o )        out_dir=$(readlink -m "$2"); shift ;;
        -t )        threads="$2"; shift ;;
        -opts )     shift; MetaSPADES_opts="$@"; break ;;
        * )         echo "Unknown option: $1"; exit 1 ;;
    esac
    shift
done

echo "Number of threads: ${threads:=4}"

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists.
validate_output_dir

# Activate conda environment
activate_env metabiome-genome-assembly

# Run metaSPAdes #

for file in "$input_dir"/*; do

    # Make sure to process only fastq, fq.gz or fastq.gz files
    if [[ "$file" == @(*_R1_*|*_1).@(fastq|fq.gz|fastq.gz) ]]; then

        echo "Performing PE assembly with files $(basename $file) and $(basename $file | forward_to_reverse)"
        file_dir="$out_dir"/$(basename "$file" | remove_forward_suffix)
        mkdir "$file_dir"
        spades.py --meta \
            -o "$file_dir" \
            -1 "$file" \
            -2 $(echo "$file" | forward_to_reverse) \
            -t "$threads" $MetaSPADES_opts # Obtain other options for metaSPAdes
    fi

done
