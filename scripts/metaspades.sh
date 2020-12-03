#!/bin/bash
# Metagenomic assembly using metaSPAdes from SPAdes-3.12.0
# Written by: Phagomica_club

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: metabiome metaspades -i <input directory> -o <output directory> [-t <threads>] [MetaSPADES_OPTIONS]"
    echo ""
    echo "Output directory will be created if it doesn't exists."
}

# Exit if command is called with no arguments
validate_arguments "$#"

while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift ;;
        -o )        out_dir=$(readlink -m "$2"); shift ;;
        -t )        threads="$2"; shift ;;
        * )         MetaSPADES_opts="$@" ;;
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
