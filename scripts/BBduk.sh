#!/bin/bash
# BBduk wrapper script for the DNAr 16S picking strategy from metagenomic samples:
# Written by: Phagomica Group
# Last updated on: 2020-08-19

set -e

function usage() {
    echo "Usage: $0 -i <input_directory> -o <output_directory> \
		-D <16S_DATABASE> [-e <conda_env>] [-t <threads>]"
    echo ""
    echo "<input_directory>  Input directory containing FASTQ files."
    echo "<out_directory> Directory in which results will be saved. This directory"
    echo "          will be created if it doesn't exists."
    echo "<conda_env>    Current conda environment."
    echo "<16S_DATABASE>    16S Database directory."
    echo ""
    echo "Options:"
    echo "-t        Number of threads to use."

}

#Saving input orders into variables:
if [[ "$#" == 0 ]]; then
    echo "No arguments given."
    usage
    exit 1
fi

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
        -e )        conda_env="$2"
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


# Output info
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=4}"
echo "Conda environment: ${conda_env:?'=conda environment not set'}"
echo "Reference database: ${database:?'=reference database not set'}"

##Matching reads against the 16S rDNA SSU from SILVA Database,

#For PE reads:

for i in "$reads_dir"/*paired_unaligned.fq.1.gz;do
  bbduk.sh in=$i in2=$(echo $i | 's/.1.gz/.2.gz/') \
	ref="$database" outm="$out_dir"/$(echo $(basename -- $i) | sed 's/.1.gz/BBduk_1_16S.fq/') \
  outm2="$out_dir"/$(echo $(basename -- $i) | sed 's/.1.gz/BBduk_2_16S.fq/') \
	outs="$out_dir"/$(echo $(basename -- $i) | sed 's/.1.gz/BBduk_single_16S.fq/') \
  stats="$out_dir"/$(echo $(basename -- $i) | sed 's/fq.1.gz/BBduk_stats.txt/') ordered=T ;done

#For SE reads:
for i in "$input_dir"/*unpaired_unaligned.fq.gz;do
  bbduk.sh in=$i ref="$database" outm="$out_dir"/$(echo $(basename -- $i) | sed 's/.fq.gz/-16S.fq/') \
	stats="$out_dir"/$(echo $(basename -- $i) | sed 's/unpaired_unaligned.fq.gz/unpaired_statistics_16S.txt/') \
  ordered=T;done

####Compressing################:
cd "$out_dir"
gzip *.fq

echo "Done."
echo "You can now use these 16S clean reads to:"
echo "- Amplicon-based analysis in QIIME2 or Mothur"
