#!/bin/bash
# Kaiju wrapper script for the taxonomic binning from metagenomic samples
# Written by: Phagomica Group
# Last updated on: 2020-08-19

set -e

function usage() {
    echo "Usage: $0 -i <input directory> -o <output directory> \
		-D <database_directory> -d <name_kaiju_database> -e <conda_env> [-t <threads>]"
    echo ""
    echo "<input_directory>  Input directory containing FASTQ files."
    echo "<output_directory> Directory in which results will be saved. This directory"
    echo "          will be created if it doesn't exists."
    echo "<conda_env>    Current conda environment."
    echo "<database_directory>    Kaiju Database directory."
    echo "<name_kaiju_database>   Kaiju Database name. Please refer to the Kaiju"
    echo "github repository in order to know the names"
    echo ""
    echo "Options:"
    echo "<threads>        Number of threads to use."
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
        -i )        input_dir=$(readlink -f "$2")
            shift
            ;;
        -o )        out_dir=$(readlink -f "$2")
            shift
            ;;
        -D )        database=$(readlink -f "$2")
            shift
            ;;
        -d )        name_kaijudb="$2"
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
echo "Kaiju reference database: ${database:?'=Kaiju database directory not set'}"
echo "Kaiju database name: ${name_kaijudb:?'=Kaiju database name not set'}"
echo "Kaiju version: $(kaiju --version)"

###----------------------Making Kaiju Database------------------------------###:
cd "$database"
kaiju-makedb -s "$name_kaijudb" -t "$threads"

###----------------------Paired-end reads-----------------------------------###:
for i in "$input_dir"/*paired_unaligned.fq.1.gz;do
  kaiju -t "$database"/*mnodes.dmp -f "$database"/*.fmi -i "$i" \
  -j $(echo $i | sed 's/.1.gz/.2.gz/') -o "$out_dir" -z "$threads";done

###----------------------Single-end reads-----------------------------------###:
for i in "$input_dir"/*unpaired_unaligned.fq.gz;do
  kaiju -t "$database"/*nodes.dmp -f "$database"/*.fmi -i "$i" \
  -o "$out_dir" -z "$threads";done
