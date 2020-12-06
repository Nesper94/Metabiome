#!/bin/bash
# HUMAnN 2.0 wrapper script
# Written by: Phagomica group
# Last updated on: 2020-08-29

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

usage() {
    echo "Usage: metabiome humann2 [Options] -i <input directory> -o <output directory>"
    echo ""
    echo "  <input directory>   directory containing .fq.gz files"
    echo "  <output directory>  directory to write results. This directory will be created if it doesn't exists."
    echo ""
    echo "Options:"
    echo "  -t  NUM         Number of threads to use"
    echo "  -h, --help      Show this help"
    echo "Note: Before running this script you should have downloaded the nucleotide and protein databases with
    humann2_databases. Run 'humann2_databases --help' for more info about this."
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
for forward in "$input_dir"/*_1_paired*.fq.gz; do
    reverse=$(echo "$forward" | sed 's/_1_paired/_2_paired/')
    TARGET="$out_dir"/tmp/$(basename "$forward" | sed 's/_1_paired/_paired-reads/')
    cat "$forward" "$reverse" > "$TARGET"
done

# For loop to run HUMAnN 2.0 on each sample
for file in "$out_dir"/tmp/*.fq.gz; do
    humann2 -i "$file" -o "$out_dir" --threads "$threads"
done

# Delete temporary directory
rm -rf "$out_dir"/tmp
