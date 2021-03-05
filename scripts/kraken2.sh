#!/bin/bash
# Kraken2 wrapper script for the taxonomic binning of reads.
# Written by: Phagomica group
# Last updated on: 2021-02-08

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat <<HELP_USAGE
Perform taxonomic classification of sequences.
Usage: metabiome kraken2 [Options] -i <input directory> -o <output directory>
Required:
  -i in_dir       Input directory containing FASTQ files."
  -o out_dir      Directory in which results will be saved. This directory"
                  will be created if it doesn't exist."

Options:
  -db database    Path to database used to assign the taxonomic labels to sequences (default: standard-kraken2-db)"
  -t  NUM         Number of threads to use (default: 1)"
  -h, --help      Show this help"
HELP_USAGE
}

# Exit if command is called with no arguments
validate_arguments "$#"

# Get input parameters
while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift ;;
        -o )        out_dir=$(readlink -m "$2"); shift ;;
        -db )       DBNAME="$2"; shift ;;
        -t )        threads="$2"; shift ;;
        * )         echo "Option '$1' not recognized" >&2; exit 1 ;;
    esac
    shift
done

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists.
validate_output_dir

# Activate conda environment
activate_env metabiome-taxonomic-binning

# Output info
echo "Input directory: $input_dir"
echo "Output directory: $out_dir"
echo "Number of threads: ${threads:=1}"
echo "Database name: ${DBNAME:=standard-kraken2-db}"
echo "Kraken2 version: $(kraken2 -v)"
echo

# Download and create database
if [[ ! -d "$DBNAME" ]]; then
    echo "Creating $DBNAME database..."
    # Create standard Kraken 2 database.
    # Installing Kraken2 through Conda using a YML file brings troubles when
    # Kraken2 tries to use rsync, for this reason we use the flag --use-ftp
    # See: https://github.com/bioconda/bioconda-recipes/issues/14076#issuecomment-474396810
    # and: https://www.biostars.org/p/423118/#428376
    kraken2-build --standard --db "$DBNAME" --use-ftp

    # Build database (yes, it seems lke we are repeating the previous command, but
    # this is a necessary step)
    echo "Building $DBNAME database..."
    kraken2-build --standard --threads "$threads" --db "$DBNAME" --use-ftp
fi

# Classification
echo "Classifying reads..."
for file in "$input_dir"/*; do
    # Make sure to process only fastq, fq.gz or fastq.gz files
    if [[ "$file" == @(*_R1*|*_1).@(fastq|fq.gz|fastq.gz) ]]; then
        forward_file="$file"
        reverse_file=$(echo "$forward_file" | forward_to_reverse)
        core_name=$(get_core_name "$forward_file")

        kraken2 --paired \
            --db "$DBNAME" \
            --threads "$threads" \
            --classified-out "$out_dir"/${core_name}_paired_classified_seqs#.fq \
            --unclassified-out "$out_dir"/${core_name}_paired_unclassified_seqs#.fq \
            --report "$out_dir"/${core_name}_paired_report.txt \
            --output "$out_dir"/${core_name}_paired_kraken2_out.tsv \
            "$forward_file" "$reverse_file"

    # Unpaired reads
    elif [[ "$file" == *unpaired*fq.gz ]]; then
        core_name=$(get_core_name "$file")

        kraken2 --db "$DBNAME" \
            --threads "$threads" \
            --classified-out "$out_dir"/${core_name}_classified_seqs.fq \
            --unclassified-out "$out_dir"/${core_name}_unclassified_seqs.fq \
            --report "$out_dir"/${core_name}_report.txt \
            --output "$out_dir"/${core_name}_kraken2_out.tsv \
            "$file"
    fi
done
