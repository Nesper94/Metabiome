#!/bin/bash
# MetaPhlAn3 wrapper script for the taxonomic profiling of metagenomic reads.
# Written by: Phagomica Group
# Last updated on: 2021-05-24

set -e
SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat <<HELP_USAGE
Profile the taxa of your metagenomic samples through MetaPhlAn3.
Usage: metabiome metaphlan3 -i <in_dir> -o <out_dir> -opts metaphlan3_options

Required:
  -i    in_dir      Input directory containing FASTQ reads files.
  -o    out_dir     Directory in which results will be saved. This directory
                    will be created if it doesn't exist.

Options:
  -d    db_dir      Directory containing the MetaPhlAn3 Database. (The database
                    will be automatically downloaded if it this directory is not set).
  -t    NUM         Number of threads to use. (default=1)
  -opts OPTIONS     MetaPhlAn3's options.
  -h, --help        Show this help.
  -hh               Show MetaPhlAn3's help message.
HELP_USAGE
}

# Exit if command is called with no arguments
validate_arguments "$#"

# Parse arguments
while (("$#")); do
    case "$1" in
        -h|--help ) usage; exit 0;;
        -hh )       activate_env metabiome-taxonomic-profiling; metaphlan -h; exit 0;;
        -i )        input_dir=$(readlink -f "$2"); shift 2;;
        -o )        out_dir=$(readlink -m "$2"); shift 2;;
        -d )        met_db=$(readlink -m "$2"); shift 2;;
        -t )        threads="$2"; shift 2;;
        -opts )     shift; metaphlan_opts="$@"; break ;;
        * )         echo "Option '$1' not recognized"; exit 1;;
    esac
done

# Verify that input directory exists
validate_input_dir

# Create output directory if it doesn't exists
validate_output_dir

# Activate conda environment
activate_env metabiome-taxonomic-profiling

# Install MetaPhlAn database
# First check if MetaPhlAn3 database directory is provided, otherwise it will be created
# and the database downloaded.
if [[ ! -d "$met_db" ]]; then
    create_dir "$out_dir" database
    met_db="$out_dir"/database
    metaphlan --install --bowtie2db "$met_db"

# Generate the index of the provided MetaPhlAn3 database
elif [[ -d "$met_db" ]]; then
    for md5_file in "$met_db"/*.md5; do
        name_idx=$(get_core_name "$md5_file")
        metaphlan --install --index "$name_idx" --bowtie2db "$met_db"
    done
fi

# Output info
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: $input_dir"
echo "Output directory: $out_dir"
echo "Number of threads: ${threads:=1}"
echo "MetaPhlAn3 database: ${met_db:?'Database not downloaded'}"
echo "MetaPhlAn3 version: $(metaphlan -v)"
echo "MetaPhlAn3 called with options: $metaphlan_opts"

# Taxonomic profiling
for file in "$input_dir"/*; do
    # Paired end reads
    if [[ "$file" == @(*_R1_*|*_1).@(fq|fastq|fq.gz|fastq.gz) ]]; then
        forward_file="$file"
        core_name=$(get_core_name "$forward_file")
        metaphlan "$forward_file",$(forward_to_reverse "$forward_file") \
            --input_type fastq  --bowtie2db "$met_db" --nproc "$threads" \
            -o "$out_dir"/$(echo "$core_name" | sed 's/_bt2/_mphlan/').txt \
            --bowtie2out "$out_dir"/$(echo "$core_name" | sed 's/_bt2/_mphlan/').sam \
            $metaphlan_opts

    # Unpaired reads
    elif [[ ! "$file" ==  *_@(*R1_*|*1.|*R2_*|*2.)* ]] && [[ "$file" == *.@(fq|fastq|fq.gz|fastq.gz) ]]; then
        unpaired_file="$file"
        core_name=$(get_core_name "$unpaired_file")
        metaphlan "$unpaired_file" --input_type fastq --bowtie2db "$met_db" \
            -o "$out_dir"/$(echo "$core_name" | sed 's/_bt2/_mphlan/').txt --nproc "$threads" \
            --bowtie2out "$out_dir"/$(echo "$core_name" | sed 's/_bt2/_mphlan/').sam \
            $metaphlan_opts

    # Files that do not match the required extension:
    elif [[ ! "$file" == *.@(fq|fastq|fq.gz|fastq.gz) ]]; then
        echo -e "$(basename -- "$file") will not be processed as is not a .fastq or .fq.gz file."
    fi
done

# Merge tables from utility scripts
cd "$out_dir"
merge_metaphlan_tables.py *.txt > merged_abundance_table.txt

echo "Done."
echo "You can now use this metagenomic profiling to:"
echo "- Play around in R or Python"
