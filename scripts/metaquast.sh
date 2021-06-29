#!/bin/bash
# metaQUAST 5.0 wrapper script
# Last updated on: 2021-05-24

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat <<HELP_USAGE
Evaluate genome assembly.
Usage: metabiome metaquast [Options] -i <in_dir> -o <out_dir> [-opts metaQUAST_OPTS]

Required:
  -i in_dir       Directory containing contigs in FASTA format
  -o out_dir      Directory to write results. This directory
                  will be created if it doesn't exists.

Options:
  -t NUM          Number of threads to use
  -opts OPTIONS   Additional options to use with metaQUAST
  -h, --help      Show this help
  -hh             Show MetaQUAST's help message.
HELP_USAGE
}

# Exit if command is called with no arguments
validate_arguments "$#"

# Parse command line arguments
while (("$#")); do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -hh )       activate_env metabiome-genome-assembly; metaquast.py -h; exit 0 ;;
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
echo "Number of threads: ${threads:=1}"

contigs=""
for file in "$input_dir"/*; do
    if [[ "$file" == *contigs*.@(fa|fasta|fna|ffn) ]]; then
        contigs="$contigs $file"
    fi
done

metaquast.py $contigs -o "$out_dir" --threads "$threads" $metaQUAST_OPTS
