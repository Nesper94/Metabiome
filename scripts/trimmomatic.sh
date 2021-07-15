#!/bin/bash
# Trimmomatic wrapper script to make quality control on FASTQ files
# Last updated on: 2021-05-24

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat <<HELP_USAGE
Perform quality trimming on Illumina sequence data.
Usage: metabiome trimmomatic [Options] -i <input directory> -o <output directory> -opts TRIMMOMATIC_OPTIONS

Required:
  -i in_dir	   Input directory containing clean FASTQ files.
  -o out_dir   	   Directory in which results will be saved. This directory
               	   will be created if it doesn't exist.

Options:
  -t NUM           Number of threads to use (default: 4)
  -opts OPTIONS    Trimmomatic's options.
  -h, --help       Show this help
  -hh              Show Trimmomatic's help message.
HELP_USAGE
}

# Exit if command is called with no arguments
validate_arguments "$#"

while (("$#")); do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -hh )       activate_env metabiome-preprocessing; trimmomatic -h; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift 2 ;;
        -o )        out_dir=$(readlink -m "$2"); shift 2 ;;
        -t )        threads="$2"; shift 2 ;;
        -opts )     shift; trimopt="$@"; break ;;
        * )         echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists.
validate_output_dir

# Activate conda environment
activate_env metabiome-preprocessing

# Output info
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: $input_dir"
echo "Output directory: $out_dir"
echo "Number of threads: ${threads:=4}"
echo "Trimmomatic version: $(trimmomatic -version)"
echo "Trimmomatic called with options: $trimopt"
echo

for file in "$input_dir"/*; do

    # Make sure to process only fastq, fq.gz or fastq.gz files
    if [[ "$file" == @(*_R1_*|*_1).@(fastq|fq.gz|fastq.gz) ]]; then
        trimmomatic PE -threads ${threads:=4} "$file" $(echo "$file" | forward_to_reverse) \
        "$out_dir"/$(basename "$file" | sed 's/\(_R1_\|_1\).*/_paired_trim_1.fq.gz/')   \
        "$out_dir"/$(basename "$file" | sed 's/\(_R1_\|_1\).*/_unpaired_trim_f.fq.gz/') \
        "$out_dir"/$(basename "$file" | sed 's/\(_R1_\|_1\).*/_paired_trim_2.fq.gz/')   \
        "$out_dir"/$(basename "$file" | sed 's/\(_R1_\|_1\).*/_unpaired_trim_r.fq.gz/') \
        $trimopt 2>&1 | tee -a "$out_dir"/trimmomatic.log

        echo # Add new line in order to make output easier to read

    elif [[ ! "$file" == *.@(fastq|fq.gz|fastq.gz) ]]; then
        echo -e "$(basename $file) will not be processed as is not a .fastq or .fq.gz file."
    fi
done

echo
echo "Done. You should now execute Bowtie2 in order to clean contaminant reads."
