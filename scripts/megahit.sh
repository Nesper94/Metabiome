#!/bin/bash
# Metagenomic assembly using MEGAHIT-1.2.9
# Written by: Phagomica_club

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat <<HELP_USAGE
Alignment of sequences derived from metagenomic samples with MEGAHIT.
Usage: metabiome megahit [Options] -i <input directory> -o <output directory> [-opts Megahit_OPTIONS]
Required:
  -i in_dir       Input directory containing FASTQ files.
  -o out_dir      Directory in which results will be saved. This directory
                  will be created if it doesn't exist.

Options:
  -t NUM          Number of threads to use (default: 4).
  -opts OPTIONS   Megahit's options.
  -h, --help      Show this help.
HELP_USAGE
}

# Exit if command is called with no arguments
validate_arguments "$#"

while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift ;;
        -o )        out_dir=$(readlink -m "$2"); shift ;;
        -t )        threads="$2"; shift ;;
        --k-list)   kmer_list="$2"; shift ;;
        -opts )     shift; MEGAHIT_opts="$@"; break ;;
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

#Output info
echo "Input directory: $input_dir"
echo "Output directory: $out_dir"
echo "Number of threads: ${threads:=4}"
echo "MEGAHIT version: $(megahit --version)"

# Run MEGAHIT on paired-end files
for file in "$input_dir"/*; do
    # Make sure to process only fastq, fq.gz or fastq.gz files
    if [[ "$file" == @(*_R1_*|*_1).@(fastq|fq.gz|fastq.gz) ]]; then

        echo "Performing PE assembly with $(basename $file) and $(basename "$file" | forward_to_reverse)"
        out_name=$(get_core_name "$file")
        megahit -o "$out_dir"/"$out_name" \
            -1 "$file" \
            -2 $(echo "$file" | forward_to_reverse) \
            -t "$threads" $MEGAHIT_opts \
            --presets meta-large # Optimization for large & complex metagenomes, like soil
    fi
done

echo "Done."
echo "You can now use assembled reads to:"
echo "- Perform taxonomic binning with kraken2.sh or MetaPhlAn3.sh"
echo "- Perform functional profiling using humann2.sh"
echo "- Extract 16S sequences with BBduk.sh"
