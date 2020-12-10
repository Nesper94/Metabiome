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
    echo "  -d db_dir                 MetaPhlan3 Database directory."
    echo ""
    echo "Options:"
    echo "  -t NUM                    Number of threads to use. (default=1)"
    echo "  -opts metaphlan_OPTIONS   metaphlan's options."
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

##----------------Installing MetaPhlAn database-----------------------------##:
##First checks if MetaPhlAn3 database already generated, otherwise it will be generated.
[ ! -f "$met_database"/mpa_v30_CHOCOPhlAn_201901.4.bt2 ] \
&& { echo "MetaPhlAn3 database needs to be downloaded" ; \
     metaphlan --install --bowtie2db "$met_database"; }

##---------------------Output info----------------------------------------##:
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=1}"
echo "MetaPhlAn3 database: ${met_database:?'Database not downloaded'}"
echo "MetaPhlAn3 version: $(metaphlan -v)"

##---------------------------MetaPhlAn profiling----------------------------##:
##----------------------PE reads--------------------------------------------##:
for foward_file in "$input_dir"/*_1_paired_;do
  metaphlan "$forward_file",$(echo "$forward_file" | sed 's/_1_/_2_/') \
  --input_type fastq  -t rel_ab_w_read_stats \
  -o "$out_dir"/$(echo $(basename -- "$forward_file") | sed 's/1_paired_bt2.fq.gz/paired_mphlan.txt/') \
  --nproc "$threads" --bowtie2out "$out_dir"/$(echo $(basename -- "$forward_file") | sed 's/1_paired_bt2.fq.gz/paired_bt2_mphlan.sam/') \
  "$metaphlan_opts"
done

##------------------------------SE reads------------------------------------##:
for unpaired_file in "$input_dir"/*_unpaired_;do
  metaphlan "$unpaired_file" --input_type fastq -t rel_ab_w_read_stats \
  -o "$out_dir"/$(echo $(basename -- "$unpaired_file") | sed 's/unpaired_bt2.fq.gz/unpaired_mphlan.txt/') \
  --nproc "$threads" --bowtie2out "$out_dir"/$(echo $(basename -- "$unpaired_file") | sed 's/unpaired_bt2.fq.gz/un_bt2_mphlan.sam/') \
  "$metaphlan_opts"
done

##-----------------------------Merging tables from utility scripts-----------##:
cd "$out_dir"
merge_metaphlan_tables.py *.txt > merged_mphlan.txt


echo "Done."
echo "You can now use this metagenomic profiling to :"
echo "- Play around in R or Python"
