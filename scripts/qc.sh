#!/bin/bash
#FastQC and MultiQC script to make a screening of FASTQ files quality
#Made by: PhagePipe 

set -e
shopt -s extglob

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: metabiome qc -i <input directory> -o <output directory> [-t <threads>] "
    echo ""
    echo "Output directory will be created if it doesn't exists."
}

# Exit if command is called with no arguments
validate_arguments "$#"

while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift ;;
        -o )        out_dir=$(readlink -f "$2"); shift ;;
        -t )        threads="$2"; shift ;;
        * )         echo "unknow option" >&2; exit 1 ;; 
    esac
    shift
done

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists.
validate_output_dir

# Activate conda environment
activate_env metabiome-preprocessing

# Output info
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: ${input_dir}"
echo "Output directory: ${out_dir}"
echo "Number of threads: ${threads:=2}"
#echo "FastQC version: $(fastqc --version)"


# Make sure to process only fastq, fq.gz or fastq.gz files

fastqc -o "$out_dir" "$input_dir"/*@(.fastq|.fastq.gz|.fq.gz)

multiqc -o "$out_dir" "$out_dir"

