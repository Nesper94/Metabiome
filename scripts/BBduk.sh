#!/bin/bash
# BBduk wrapper script for the DNAr 16S picking strategy from metagenomic samples:
# Written by: Phagomica Group
# Last updated on: 2020-27-10

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: $0 -i <input_directory> -o <output_directory> -D <16S_DATABASE> [-t <threads>]"
    echo ""
    echo "Options:"
    echo "<input_directory>  Input directory containing FASTQ files."
    echo "<output_directory> Directory in which results will be saved. This directory"
    echo "will be created if it doesn't exist."
    echo "<16S_DATABASE> 16S Database directory."
    echo "<threads>  Number of threads to use. (optional)"

}

# Exit if command is called with no arguments
validate_arguments "$#"

##---------------------Save input parameters into variables------------------##:
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
        -D )        database=$(readlink -f "$2")
            shift
            ;;
        -t )        threads="$2"
           shift
            ;;
        * )        echo "Option '$1' not recognized"; exit 1
            ;;
    esac
    shift
done

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists.
validate_output_dir

##----------------------------Output info------------------------------------##:
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=4}"
echo "Reference database: ${database:?'=reference database not set'}"

##-----------------------Activate conda environment--------------------------##:
activate_env picking16S

##------------Match reads against the 16S rDNA SSU from SILVA Database-------##:

for i in "$input_dir"/*1_paired_bt2.fq.gz; do
    echo $i
    bbduk.sh in=$i in2=$(echo $i | 's/_1_paired_bt2/_2_paired_bt2/') \
	ref="$database" outm="$out_dir"/$(echo $(basename -- $i) | sed 's/1_paired_bt2.fq.gz/1_paired_bbdk.fq/') \
    outm2="$out_dir"/$(echo $(basename -- $i) | sed 's/1_paired_bt2.fq.gz/2_paired_bbdk.fq/') \
	outs="$out_dir"/$(echo $(basename -- $i) | sed 's/1_paired_bt2.fq.gz/singletons_bbdk.fq/') \
    stats="$out_dir"/$(echo $(basename -- $i) | sed 's/1_paired_bt2.fq.gz/pe_bbdk_summary.txt/') ordered=T
done

##----------------------------For SE reads-----------------------------------##:

for i in "$input_dir"/*unpaired_bt2.fq.gz; do
    bbduk.sh in=$i ref="$database" \
    outm="$out_dir"/$(echo $(basename -- $i) | sed 's/unpaired_bt2.fq.gz/unpaired_bbdk.fq/') \
	stats="$out_dir"/$(echo $(basename -- $i) | sed 's/unpaired_bt2.fq.gz/unpaired_bbdk_summary.txt/') \
    ordered=T
done

##-------------------------------Compress output-----------------------------##:
cd "$out_dir"
gzip *.fq

echo "Done."
echo "You can now use these 16S clean reads to:"
echo "- Amplicon-based analysis in QIIME2 or Mothur"
