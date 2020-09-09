#!/bin/bash
# MetaPhlAn3 wrapper script for the taxonomic profiling of reads.
# Written by: Phagomica Group
# Last updated on: 2020-08-28

##------------------------------Checking the input---------------------------##:
set -e
function usage() {
    echo "Usage: $0 -i <input directory> -o <output directory> -e <conda_env> \
    -d <database_folder> [-t <threads>] "
    echo "Output directory will be created if it doesn't exists."
}

if [[ "$#" == 0 ]]; then
    echo "No arguments given."
    usage
    exit 1
fi

while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0
            ;;
        -i )        reads_dir=$(readlink -f "$2")
            shift
            ;;
        -o )        out_dir=$(readlink -f "$2")
            shift
            ;;
        -e )        conda_env="$2"
            shift
            ;;
        -d )        met_database=$(readlink -f "$2")
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

# Verify that input directory exists
if [ ! -d "$input_dir" ]; then
   echo "$0: Error: $input_dir is not a valid directory."
   exit 1
fi
# Create output directory if it doesn't exists.
if [[ ! -d "$out_dir" ]]; then
    mkdir "$out_dir"
fi

##----------------Installing MetaPhlAn database--------------------------##:

metaphlan --install --bowtie2db "$met_database"

# Output info
echo "Input directory: ${reads_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=4}"
echo "MetaPhlAn3 conda environment: ${conda_env:?'conda environment not set'}"
echo "MetaPhlAn3 database: ${met_database:?'Database not downloaded'}"
echo "MetaPhlAn3 version: $(metaphlan -v)"

##---------------------------MetaPhlAn profiling----------------------------##:
##----------------------PE reads-----------------------------------##:
for i in "$reads_dir"/*paired_unaligned.fq.1.gz;do
  metaphlan $i,$(echo $i | sed 's/.1.gz/.2.gz/') \
  --input_type fastq --add_viruses  -t rel_ab_w_read_stats --unknown_estimation \
  -o "$out_dir"/$(echo $(basename -- $i) | sed 's/.fq.1.gz/metaphlan_out.txt/') \
  --nproc "$threads" --bowtie2out "$out_dir"/$(echo $(basename -- $i) | sed 's/.fq.1.gz/bowtie2_out/') \
  ;done

##---------------------SE reads------------------------------------##:
for i in "$reads_dir"/*unpaired_unaligned.fq.gz;do
  metaphlan $i --input_type fastq --add_viruses --unknown_estimation \
  -t rel_ab_w_read_stats -o "$out_dir"/$(echo $(basename -- $i) | sed 's/.fq.gz/metaphlan_out.txt/') \
  --nproc "$threads" --bowtie2out "$out_dir"/$(echo $(basename -- $i) | sed 's/.fq.gz/bowtie2_out/') ;done
