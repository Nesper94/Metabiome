#!/bin/bash
# MetaPhlAn3 wrapper script for the taxonomic profiling of reads.
# Written by: Phagomica Group
# Last updated on: 2020-10-28

##------------------------------Checking the input---------------------------##:
set -e
SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: metabiome metaphlan3 -i <input_directory> -o <output_directory> -d <database_directory> [-t <threads>] [metaphlan_OPTIONS]"
    echo ""
    echo "Options:"
    echo "<input_directory>  Input directory containing clean FASTQ files."
    echo "<output_directory> Directory in which results will be saved. This directory"
    echo "will be created if it doesn't exist."
    echo "<database_directory>  MetaPhlan3 Database directory."
    echo "<threads>  Number of threads to use. (default=1)"
    echo "<metaphlan_OPTIONS> metaphlan's options (optional). Make sure to enclose metaphlan_OPTIONS within quotation marks"
}


##--------------------------Exiting if input files are missing---------------##:
validate_arguments "$#"


##----------------------Saving input orders into variables------------------##:
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
        -d )        met_database=$(readlink -f "$2")
            shift
            ;;
        -t )        threads="$2"
           shift
            ;;
        * )         metaphlan_opts="$@"
            ;;
        * )        echo "Option '$1' not recognized"; exit 1
            ;;
    esac
    shift
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
for i in "$input_dir"/*_1_paired_bt2.fq.gz;do
  metaphlan $i,$(echo $i | sed 's/_1_/_2_/') \
  --input_type fastq  -t rel_ab_w_read_stats --unknown_estimation \
  -o "$out_dir"/$(echo $(basename -- $i) | sed 's/1_paired_bt2.fq.gz/paired_mphlan.txt/') \
  --nproc "$threads" --bowtie2out "$out_dir"/$(echo $(basename -- $i) | sed 's/1_paired_bt2.fq.gz/paired_bt2mphlan.sam/') \
  "${metaphlan_opts:=''}"
done

##------------------------------SE reads------------------------------------##:
for i in "$input_dir"/*unpaired_bt2.fq.gz;do
  metaphlan $i --input_type fastq -t rel_ab_w_read_stats \
  -o "$out_dir"/$(echo $(basename -- $i) | sed 's/unpaired_bt2.fq.gz/unpaired_mphlan.txt/') \
  --nproc "$threads" --bowtie2out "$out_dir"/$(echo $(basename -- $i) | sed 's/unpaired_bt2.fq.gz/un_bt2mphlan.sam/') \
  "${metaphlan_opts:=''}"
done

##-----------------------------Merging tables from utility scripts-----------##:

cd $out_dir
merge_metaphlan_tables.py *.txt > merged_mphlan.txt


echo "Done."
echo "You can now use this metagenomic profiling to :"
echo "- Play around in R or Python"
