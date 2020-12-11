#!/bin/bash
# MetaPhlAn3 wrapper script for the taxonomic profiling of reads.
# Written by: Phagomica Group
# Last updated on: 2020-10-28

##------------------------------Checking the input---------------------------##:
set -e
SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: metabiome metaphlan3 -i <in_dir> -o <out_dir> -d <db_dir>"
    echo ""
    echo "Mandatory:"
    echo "  -i in_dir                 Input directory containing clean FASTQ files."
    echo "  -o out_dir                Directory in which results will be saved. This directory"
    echo "                            will be created if it doesn't exist."
    echo "  -d db_dir                 MetaPhlAn3 Database directory."
    echo ""
    echo "Options:"
    echo "  -t NUM                    Number of threads to use. (default=1)"
    echo "  -opts MetaPhlAn_OPTIONS   MetaPhlAn3's options."
    echo "  -h, --help                Show this help"
}

##--------------------------Exiting if input files are missing---------------##:
validate_arguments "$#"


##----------------------Saving input orders into variables------------------##:
while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0;;
        -i )        input_dir=$(readlink -f "$2"); shift 2;;
        -o )        out_dir=$(readlink -m "$2");shift 2;;
        -d )        met_database=$(readlink -f "$2"); shift 2;;
        -t )        threads="$2"; shift 2;;
        -opts )     shift; metaphlan_opts="$@"; break ;;
        * )         echo "Option '$1' not recognized"; exit 1;;
    esac
done

##----------------Verify that input directory exists------------------------##:
validate_input_dir
##---------------Create output directory if it doesn't exists---------------##:
validate_output_dir
##----------------Activate conda environment--------------------------------##:
activate_env metabiome-metaphlan

##--------------------Install MetaPhlAn database-----------------------------##:
# First checks if MetaPhlAn3 database already generated, otherwise it will be generated.
[ ! -f "$met_database"/mpa_v30_CHOCOPhlAn_201901.4.bt2 ] \
&& { echo "MetaPhlAn3 database needs to be downloaded" ; \
     metaphlan --install --bowtie2db "$met_database"; }

##---------------------Output info----------------------------------------##:
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: $input_dir"
echo "Output directory: $out_dir"
echo "Number of threads: ${threads:=1}"
echo "MetaPhlAn3 database: ${met_database:?'Database not downloaded'}"
echo "MetaPhlAn3 version: $(metaphlan -v)"

##---------------------------MetaPhlAn profiling----------------------------##:
for file in "$input_dir"/*; do
##-----------------------------Paired End reads------------------------------##:
    if [[ "$file" == @(*_R1_*|*_1).@(fastq|fq.gz|fastq.gz) ]]; then
        forward_file="$file"
        metaphlan "$forward_file",$(echo "$forward_file" | forward_to_reverse) \
            --input_type fastq -t rel_ab_w_read_stats \
            -o "$out_dir"/$(get_core_name "$forward_file" | sed 's/_bt2//')_mphlan.txt \
            --nproc "$threads" \
            --bowtie2out "$out_dir"/$(get_core_name "$forward_file" | sed 's/_bt2//')_mphlan.sam \
            "$metaphlan_opts"
##------------------------------Single End reads-----------------------------##:
    elif [[ "$file" == *_unpaired_* ]]; then
        unpaired_file="$file"
        metaphlan "$unpaired_file" --input_type fastq -t rel_ab_w_read_stats \
            -o "$out_dir"/$(get_core_name "$unpaired_file" | sed 's/_bt2//')_mphlan.txt \
            --nproc "$threads" \
            --bowtie2out "$out_dir"/$(get_core_name "$unpaired_file" | sed 's/_bt2//')_mphlan.sam \
            "$metaphlan_opts"
    fi
done

##-----------------------Merge tables from utility scripts-------------------##:
cd "$out_dir"
merge_metaphlan_tables.py *.txt > merged_mphlan.txt

echo "Done."
echo "You can now use this metagenomic profiling to:"
echo "- Play around in R or Python"
