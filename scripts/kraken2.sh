#!/bin/bash
# Kraken2 wrapper script for the taxonomic binning of reads.
# Written by: Juan C. Arboleda R.
# Last updated on: 2020-08-19

set -e

function usage() {
    echo "Usage: $0 -i <input directory> -o <output directory> [-db <database name>] [-t <threads>]"
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
        -db )       DBNAME="$2"
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

# Output info
echo "Input directory: ${reads_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=4}"
echo "Database name: ${DBNAME:=standard-kraken2-db}"
echo "Kraken2 version: $(kraken2 -v)"

if [[ ! -d "$out_dir" ]]; then  # Create output directory if it doesn't exists.
    mkdir "$out_dir"
fi

if [[ "$DBNAME" == "standard-kraken2-db" ]]; then
    echo "Creating standard Kraken2 database..."
    kraken2-build --standard --db "$DBNAME" # Create standard Kraken 2 database. If -db is not set, then set as DBNAME=SILVA

    # Build database
    echo "Building standard Kraken2 database..."
    kraken2-build --standard --threads "$threads" --db $DBNAME
fi

# Classification
cd "$reads_dir"
# Paired reads
echo "Classifying paired reads..."
for file in *f-paired*.fq.gz; do
    kraken2 --paired --db "$DBNAME" --threads "$threads" --classified-out "$out_dir"/${file}-classified-seqs#.fq \
    --unclassified-out "$out_dir"/${file}-unclassified-seqs#.fq \
    --report "$out_dir"/${file}-report.txt \
    "$file" $(echo "$file" | sed 's/f-paired/r-paired/')
done

# Unpaired reads
echo "Classifying unpaired reads..."
for file in *unpaired.fq.gz; do
    kraken2 --db "$DBNAME" --threads "$threads" --classified-out "$out_dir"/${file}-classified-seqs.fq \
    --unclassified-out "$out_dir"/${file}-unclassified-seqs.fq \
    --report "$out_dir"/${file}-report.txt \
    "$file"
done
