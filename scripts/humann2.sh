#!/bin/bash
# HUMAnN 2.0 wrapper script
# Written by: Phagomica group
# Last updated on: 2020-08-29

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

usage() {
    echo "Usage: metabiome humann2 [Options] -i <input directory> -o <output directory>"
    echo
    echo "Required:"
    echo "  -i in_dir   Directory containing .fq.gz files"
    echo "  -o out_dir  Directory to write results. This directory"
    echo "              will be created if it doesn't exists."
    echo
    echo "Options:"
    echo "  -t  NUM     Number of threads to use"
    echo "  -h, --help  Show this help"
    echo
    echo "Note: Before running this script you should have downloaded the nucleotide and"
    echo "protein databases with humann2_databases. Run 'humann2_databases --help' for"
    echo "more info about this."
}

# Exit if command is called with no arguments
validate_arguments "$#"

# Parse command line arguments
while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift ;;
        -o )        out_dir=$(readlink -m "$2"); shift ;;
        -t )        threads="$2"; shift ;;
        * )         echo "Option '$1' not recognized"; exit 1 ;;
    esac
    shift
done

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists.
validate_output_dir

# Activate conda environment
activate_env metabiome-humann2

# Useful output info
echo "Input directory: ${input_dir}"
echo "Output directory: ${out_dir}"
echo "Number of threads: ${threads:=1}"
echo "HUMAnN 2.0 version: $(humann2 --version)"

# Create temporary directory to save paired reads
mkdir "$out_dir"/tmp

# Cat paired-end reads and save to tmp/.
# HUMAnN2 doesn't use paired-end information: https://forum.biobakery.org/t/paired-end-files-humann2/396/4
for file in "$input_dir"/*; do

    # Make sure to process only fastq, fq.gz or fastq.gz files
    if [[ "$file" == @(*_R1*|*_1).@(fastq|fq.gz|fastq.gz) ]]; then
        forward_file="$file"
        reverse=$(echo "$forward_file" | forward_to_reverse)
        TARGET="$out_dir"/tmp/$(basename "$forward_file" | remove_forward_suffix)
        cat "$forward_file" "$reverse" > "$TARGET"

    elif [[ ! "$file" == *.@(fastq|fq.gz|fastq.gz) ]]; then
        echo -e "$(basename $file) will not be processed as is not a .fastq or .fq.gz file."
    fi
done

# For loop to run HUMAnN 2.0 on each sample
for file in "$out_dir"/tmp/*; do
    if [[ "$file" == *.@(fastq|fq.gz|fastq.gz) ]]; then
        humann2 -i "$file" -o "$out_dir" --threads "$threads"
    fi
done

# Delete temporary directory
rm -rf "$out_dir"/tmp
