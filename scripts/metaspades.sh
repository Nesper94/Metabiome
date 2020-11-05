#!/bin/bash
# Metagenomic assembly using metaSPAdes from SPAdes-3.12.0
# Written by: Phagomica_club

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: $0 -i <input directory> -o <output directory> [-t <threads>] [MetaSPADES_OPTIONS]"
    echo ""
    echo "WARNING: Make sure to enclose MetaSPADES_OPTIONS within quotation marks."
    echo "Output directory will be created if it doesn't exists."
}

# Exit if command is called with no arguments
validate_arguments "$#"

while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0
            ;;
        -i )        input_dir=$(readlink -f "$2")
            shift
            ;;
        -o )        out_dir=$(readlink -f "$2")
            shift
            ;;
        -t )        threads="$2"
            shift
            ;;
        * )         MetaSPADES_opts="$@"
            ;;
    esac
    shift
done

echo "Number of threads: ${threads:=4}"

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists.
validate_output_dir

# Activate conda environment
activate_env assembly

# Run metaSPAdes #

forward_file_suffix=1_paired_bt2.fq.gz
reverse_file_suffix=2_paired_bt2.fq.gz

for forward_file in "$input_dir"/*$forward_file_suffix; do

	echo "Performing PE assembly with files $(basename $forward_file) and $(basename $forward_file | sed 's/_1_bt2.fq.gz/_2_bt2.fq.gz/')"
	spades.py --meta \
    -o "$out_dir" \
    -1 "$forward_file" \
    -2 $(echo "$forward_file" | sed "s/$forward_file_suffix/$reverse_file_suffix/") `# Reverse sequences` \
    -t "$threads" #"${MetaSPADES_opts:=''}" # Obtain other options for metaSPAdes

done
