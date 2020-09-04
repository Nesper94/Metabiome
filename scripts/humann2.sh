#!/bin/bash
# HUMAnN 2.0 wrapper script
# Written by: Phagomica group
# Last updated on: 2020-08-29

set -e

function usage() {
    echo "Usage: $0 -i <input directory> -o <output directory> [-t <threads>]"
    echo ""
    echo "  <input directory>   directory containing .fq.gz files"
    echo "  <output directory>  directory to write results. This directory will be created if it doesn't exists."
    echo ""
    echo "Options:"
    echo "  -t              number of threads to use"
    echo "  -h, --help      show this help"
    echo "Note: Before running this script you should have downloaded the nucleotide and protein databases with
    humann2_databases. Run 'humann2_databases --help' for more info about this."
}

if [[ "$#" == 0 ]]; then
    echo "No arguments given."
    usage
    exit 1
fi

# Parse command line arguments
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
        -t )        threads="$2"
            shift
            ;;
        * )         echo "Option '$1' not recognized"; exit 1
            ;;
    esac
    shift
done

# Useful output info
echo "Input directory: ${reads_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=4}"
echo "HUMAnN 2.0 version: $(humann2 --version)"

# Verify that reads directory exists
if [ ! -d "$reads_dir" ]; then
   echo "$0: Error: $reads_dir is not a valid directory."
   exit 1
fi

# Cat paired-end reads and delete original files.
# HUMAnN2 doesn't use paired-end information: https://forum.biobakery.org/t/paired-end-files-humann2/396/4
for forward in "$reads_dir"/*f-paired*.fq.gz; do
    reverse=$(echo "$forward" | sed 's/f-paired/r-paired/')
    cat "$forward" "$reverse" > $(echo "$forward" | sed 's/f-paired/paired-reads/') && rm "$forward" "$reverse"
done

# For loop to run HUMAnN 2.0 on each sample
for file in "$reads_dir"/*.fq.gz; do
    humann2 -i "$file" -o "$out_dir" --threads "$threads"
done
