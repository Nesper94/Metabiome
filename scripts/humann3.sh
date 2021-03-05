#!/bin/bash
# HUMAnN 3.0 wrapper script
# Written by: Estefany Lorenzana, Cristian Grisales and Juan Camilo Arboleda
# Last updated on: 2021-02-15

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

usage() {
cat <<EOF
Profile the abundance of microbial metabolic pathways and molecular functions.
Usage: metabiome humann [Options] -i <in_dir> -o <out_dir>
Required:
  -i in_dir   Directory containing FASTQ files
  -o out_dir  Directory to write results. This directory
              will be created if it doesn't exists.

Options:
  -t  NUM     Number of threads to use
  -h, --help  Show this help

Note: Before running this script you should have downloaded the nucleotide and
protein databases with humann_databases. Activate environment with
'conda activate metabiome-taxonomic-profiling' and then run
'humann_databases --help' for more info about this.
EOF
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
activate_env metabiome-humann

# Useful output info
echo "Input directory: ${input_dir}"
echo "Output directory: ${out_dir}"
echo "Number of threads: ${threads:=1}"
echo "HUMAnN version: $(humann --version)"

# Create temporary directory to save paired reads
mkdir "$out_dir"/tmp

# Cat paired-end reads and save to tmp/.
# HUMAnN doesn't use paired-end information, however it is useful to convert the
# paired reads to a single input file: https://forum.biobakery.org/t/paired-end-files-humann2/396/4
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

# For loop to run HUMAnN 3.0 on each sample
for file in "$out_dir"/tmp/*; do
    if [[ "$file" == *.@(fastq|fq.gz|fastq.gz) ]]; then
        humann --input "$file" --output "$out_dir" --threads "$threads"
    fi
done

# Delete temporary directory
rm -rf "$out_dir"/tmp
