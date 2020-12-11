#!/bin/bash
# BBduk wrapper script for the DNAr 16S picking strategy from metagenomic samples:
# Written by: Phagomica Group
# Last updated on: 2020-12-11

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: metabiome bbduk -i <in_dir> -o <out_dir> -D <16S_db> [BBduk_OPTIONS]"
    echo
    echo "Required:"
    echo "  -i in_dir             Input directory containing clean FASTQ files."
    echo "  -o out_dir            Directory in which results will be saved. This directory"
    echo "                        will be created if it doesn't exist."
    echo "  -D 16S_db             16S Database directory."
    echo
    echo "Options:"
    echo "  -t NUM                Number of threads to use. (default=1)"
    echo "  -opts BBduk_OPTIONS   BBduk's options."
    echo "  -h, --help            Show this help"
}

# Exit if command is called with no arguments
validate_arguments "$#"

while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift 2 ;;
        -o )        out_dir=$(readlink -m "$2"); shift 2 ;;
        -D )        database=$(readlink -f "$2"); shift 2 ;;
        -t )        threads="$2"; shift 2 ;;
        -opts )     shift;bbduk_opts="$@"; break ;;
        * )         echo "Option '$1' not recognized"; exit 1 ;;
    esac
done

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists.
validate_output_dir

##----------------------------Output info------------------------------------##:
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: $input_dir"
echo "Output directory: $out_dir"
echo "Number of threads: ${threads:=1}"
echo "Reference database: ${database:?'=reference database not set'}"
echo "BBduk called with options: $bbduk_opts"

##-----------------------Activate conda environment--------------------------##:
activate_env metabiome-picking16S

##------------Match reads against the 16S rDNA SSU from SILVA Database-------##:
for file in "$input_dir"/*1_paired*; do
    # Paired end reads
    if [[ "$file" == @(*_R1_*|*_1).@(fastq|fq.gz|fastq.gz) ]]; then
        forward_file="$file"
        core_name=$(get_core_name "$forward_file")
        bbduk.sh in="$forward_file" in2= $(echo "$forward_file" | forward_to_reverse) \
            ref="$database" \
            outm="$out_dir"/$(basename -- "$forward_file" | sed 's/_bt2/_bbdk/') \
            outm2="$out_dir"/$(echo "$core_name" | sed 's/_bt2//')_bbdk_2.fq \
            outs="$out_dir"/$(echo "$core_name" | sed 's/_bt2// ; s/_paired//')_singletons_bbdk.fq \
            stats="$out_dir"/$(echo "$core_name" | sed 's/_bt2//')_paired_bbdk_summary.txt \
            "$bbduk_opts"
    # Single end reads
    elif [[ "$file" == *_unpaired_* ]]; then
        unpaired_file="$file"
        core_name=$(get_core_name "$unpaired_file")
        bbduk.sh in="$unpaired_file" ref="$database" \
            outm="$out_dir"/$(echo "$core_name" | sed 's/_bt2//')_bbdk.fq \
            stats="$out_dir"/$(echo "$core_name" | sed 's/_bt2//')_bbdk_summary.txt \
            "$bbduk_opts"
    fi
done

##-------------------------------Compress output-----------------------------##:
cd "$out_dir"
gzip *.fq

echo "Done."
echo "You can now use these 16S clean reads to:"
echo "- Perform amplicon-based analysis in Qiime2 or Mothur"
