#!/bin/bash
# Metagenomic assembly using metaSPAdes from SPAdes-3.12.0
# Written by: Phagomica_club
# Last updated on: 2021-05-24

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat <<HELP_USAGE
Assembly of sequences derived from metagenomic samples with metaSPAdes.
Usage: metabiome metaspades [Options] -i <input directory> -o <output directory> [-opts MetaSPADES_OPTIONS]

Required:
  -i in_dir       Input directory containing FASTQ files.
  -o out_dir      Directory in which results will be saved. This directory
                  will be created if it doesn't exist.

Options:
  -t NUM          Number of threads to use (default: 4).
  -opts OPTIONS   MetaSPADES's options.
  -h, --help      Show this help.
  -hh             Show metaSPAdes's help message.
HELP_USAGE
}

# Exit if command is called with no arguments
validate_arguments "$#"

while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -hh )       activate_env metabiome-genome-assembly; metaspades.py -h; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift ;;
        -o )        out_dir=$(readlink -m "$2"); shift ;;
        -t )        threads="$2"; shift ;;
        -opts )     shift; MetaSPADES_opts="$@"; break ;;
        * )         echo "Unknown option: $1"; exit 1 ;;
    esac
    shift
done

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists.
validate_output_dir

# Activate conda environment
activate_env metabiome-genome-assembly

# Output info
echo "Input directory: $input_dir"
echo "Output directory: $out_dir"
echo "Number of threads: ${threads:=4}"
echo "SPAdes version: $(spades.py --version)"

# Run metaSPAdes
for file in "$input_dir"/*; do
    # Make sure to process only fastq, fq.gz or fastq.gz files
    if [[ "$file" == @(*_R1_*|*_1).@(fastq|fq.gz|fastq.gz) ]]; then

        echo "Performing PE assembly with files $(basename $file) and $(basename $file | forward_to_reverse)"
        file_dir="$out_dir"/$(get_core_name "$file")
        mkdir "$file_dir"
        spades.py --meta \
            -o "$file_dir" \
            -1 "$file" \
            -2 $(echo "$file" | forward_to_reverse) \
            -t "$threads" $MetaSPADES_opts # Obtain other options for metaSPAdes
    fi
done

echo "Done."
echo "You can now use assembled reads to:"
echo "- Perform taxonomic binning with metabiome kraken2 or metabiome metaphlan3"
echo "- Perform functional profiling using metabiome humann"
echo "- Extract 16S sequences with metabiome bbduk"
