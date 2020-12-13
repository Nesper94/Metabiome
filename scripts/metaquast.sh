#!/bin/bash
# metaQUAST 5.0 wrapper script
# Last updated on: 2020-12-13

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

usage() {
    echo "Usage: metabiome metaquast [Options] -i <in_dir> -o <out_dir> [-opts metaQUAST_OPTS]"
    echo
    echo "Required:"
    echo "  -i in_dir       Directory containing contigs in FASTA format"
    echo "  -o out_dir      Directory to write results. This directory"
    echo "                  will be created if it doesn't exists."
    echo
    echo "Options:"
    echo "  -t NUM          Number of threads to use"
    echo "  -opts OPTIONS   Additional options to use with metaQUAST"
    echo "  -h, --help      Show this help"
}

# Exit if command is called with no arguments
validate_arguments "$#"

# Parse command line arguments
while (("$#")); do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift 2;;
        -o )        out_dir=$(readlink -m "$2"); shift 2;;
        -t )        threads="$2"; shift 2;;
        -opts )     shift; metaQUAST_OPTS="$@"; break ;;
        * )         echo "Option '$1' not recognized"; exit 1 ;;
    esac
done

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists.
validate_output_dir

# Activate conda environment
activate_env metabiome-genome-assembly

# Useful output info
echo "Input directory: $input_dir"
echo "Output directory: $out_dir"
echo "Number of threads: ${threads:=4}"

contigs=""
for file in "$input_dir"/*; do
    contigs="$contigs $file"
done

metaquast.py -o "$out_dir" $metaQUAST_OPTS "$contigs"
