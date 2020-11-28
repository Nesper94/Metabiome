#!/bin/bash
# Metagenomic assembly using MEGAHIT-1.2.9
# Written by: Phagomica_club

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: metabiome megahit -i <input directory> -o <output directory> [-t <threads>] [OPTIONS]"
    echo
    echo "Output directory will be created if it doesn't exists."
}

# Exit if command is called with no arguments
validate_arguments "$#"

# Activate conda environment
activate_env metabiome-genome-assembly

while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift ;;
        -o )        out_dir=$(readlink -f "$2"); shift ;;
        -t )        threads="$2"; shift ;;
        --k-list)   kmer_list="$2"; shift ;;
        * )         MEGAHIT_opts="$@" ;;
    esac
    shift
done

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists.
validate_output_dir

#Output info
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=1}"
echo "MEGAHIT version: $(megahit --version)"

# Run MEGAHIT on paired-end files

for file in "$input_dir"/*; do

    # Make sure to process only fastq, fq.gz or fastq.gz files
    if [[ "$file" == @(*_R1_*|*_1).@(fastq|fq.gz|fastq.gz) ]]; then

        echo "Performing PE assembly with $(basename $file) and $(basename "$file" | forward_to_reverse)"
        out_name=$(basename "$file" | remove_forward_suffix)
        megahit -o "$out_dir"/"$out_name" \
            -1 "$file" \
            -2 $(echo "$file" | forward_to_reverse) \
            -t "$threads" $MEGAHIT_opts \
            --presets meta-large # Optimization for large & complex metagenomes, like soil
    fi

done
