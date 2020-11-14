#!/bin/bash
# Kraken2 wrapper script for the taxonomic binning of reads.
# Written by: Phagomica group
# Last updated on: 2020-08-19

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: $0 [Options] -i <input directory> -o <output directory>"
    echo ""
    echo "Options:"
    echo "  -db database    Path to database used to assign the taxonomic labels to sequences (default: standard-kraken2-db)"
    echo "  -t  NUM         Number of threads to use (default: 1)"
    echo "  -h, --help      Show this help"
    echo ""
    echo "Output directory will be created if it doesn't exists."
}

# Exit if command is called with no arguments
validate_arguments "$#"

# Get input parameters
while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift ;;
        -o )        out_dir=$(readlink -f "$2"); shift ;;
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
activate_env binning

# Output info
echo "Input directory: ${input_dir}"
echo "Output directory: ${out_dir}"
echo "Number of threads: ${threads:=1}"
echo "Database name: ${DBNAME:=standard-kraken2-db}"
echo "Kraken2 version: $(kraken2 -v)"

# Create standard database
if [[ "$DBNAME" == "standard-kraken2-db" ]] && [[ ! -d "$DBNAME" ]]; then
    echo "Creating standard Kraken2 database..."
    # Create standard Kraken 2 database.
    # Installing Kraken2 through Conda using a YML file brings troubles when
    # Kraken2 tries to use rsync, for this reason we use the flag --use-ftp
    # See: https://github.com/bioconda/bioconda-recipes/issues/14076#issuecomment-474396810
    # and: https://www.biostars.org/p/423118/#428376
    kraken2-build --standard --threads "$threads" --db "$DBNAME" --use-ftp
fi

# Classification

FORWARD_FILE_SUFFIX=1_paired_bt2.fq.gz
REVERSE_FILE_SUFFIX=2_paired_bt2.fq.gz

# Paired reads
echo "Classifying paired reads..."
for forward_file in "$input_dir"/*"$FORWARD_FILE_SUFFIX"; do
    kraken2 --paired \
    --db "$DBNAME" \
    --threads "$threads" \
    --classified-out "$out_dir"/${forward_file}-classified-seqs#.fq \
    --unclassified-out "$out_dir"/${forward_file}-unclassified-seqs#.fq \
    --report "$out_dir"/${forward_file}-report.txt \
    "$forward_file" $(echo "$forward_file" | sed "s/$FORWARD_FILE_SUFFIX/$REVERSE_FILE_SUFFIX/")
done

# Unpaired reads
echo "Classifying unpaired reads..."
for file in "$input_dir"/*unpaired*fq.gz; do
    kraken2 --db "$DBNAME" \
    --threads "$threads" \
    --classified-out "$out_dir"/${file}-classified-seqs.fq \
    --unclassified-out "$out_dir"/${file}-unclassified-seqs.fq \
    --report "$out_dir"/${file}-report.txt \
    "$file"
done
