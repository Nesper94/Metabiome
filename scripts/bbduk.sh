#!/bin/bash
# BBDuk wrapper script to extract 16S rDNA sequences from metagenomic samples.
# Written by: Phagomica Group
# Last updated on: 2021-05-24

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat <<HELP_USAGE
Extract 16S rDNA sequences from metagenomic samples with BBDuk.
Usage: metabiome bbduk [options] -i <in_dir> -o <out_dir> -D <16S_db> -opts bbduk_opts

Required:
  -i in_dir         Input directory containing FASTQ files.
  -o out_dir        Directory in which results will be saved. This directory will be
                    created if it doesn't exist.
  -D 16S_db         16S rRNA reference database in fasta format.

Options:
  -t NUM            Number of threads to use. (default=1)
  -opts OPTIONS     BBDuk's options.
  -h, --help        Show this help.
  -hh               Show BBDuk's help message.
HELP_USAGE
}

# Exit if command is called with no arguments
validate_arguments "$#"

# Parse command line arguments
while (("$#")); do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -hh )       activate_env metabiome-picking16S; bbduk.sh -h; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift 2 ;;
        -o )        out_dir=$(readlink -m "$2"); shift 2 ;;
        -D )        database=$(readlink -m "$2"); shift 2 ;;
        -t )        threads="$2"; shift 2 ;;
        -opts )     shift; bbduk_opts="$@"; break ;;
        * )         echo "Option '$1' not recognized"; exit 1 ;;
    esac
done

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists
validate_output_dir

# Output info
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: $input_dir"
echo "Output directory: $out_dir"
echo "Number of threads: ${threads:=1}"
echo "Reference database: ${database:?'reference database not set'}"
echo "BBDuk called with options: $bbduk_opts"

# Activate conda environment
activate_env metabiome-extract16S

# Match reads against the 16S rDNA SSU from SILVA Database
for file in "$input_dir"/*; do
    # Paired end reads
    if [[ "$file" == @(*_R1_*|*_1).@(fq|fastq|fq.gz|fastq.gz) ]]; then
        forward_file="$file"
        core_name=$(get_core_name "$forward_file")
        bbduk.sh in1="$forward_file" in2=$(echo "$forward_file" | forward_to_reverse) \
            ref="$database" \
            outm="$out_dir"/$(basename -- "$forward_file" | sed 's/_bt2/_bbdk/') \
            outm2="$out_dir"/$(echo "$core_name" | sed 's/_bt2//')_bbdk_2.fq \
            outs="$out_dir"/$(echo "$core_name" | sed 's/_bt2// ; s/_paired//')_singletons_bbdk.fq \
            stats="$out_dir"/$(echo "$core_name" | sed 's/_bt2//')_bbdk_summary.txt \
            $bbduk_opts

    # Unpaired reads
    elif [[ ! "$file" ==  *_@(*R1_*|*1.|*R2_*|*2.)* ]] && [[ "$file" == *.@(fq|fastq|fq.gz|fastq.gz) ]]; then
        unpaired_file="$file"
        core_name=$(get_core_name "$unpaired_file")
        bbduk.sh in="$unpaired_file" ref="$database" \
            outm="$out_dir"/$(echo "$core_name" | sed 's/_bt2//')_bbdk.fq \
            stats="$out_dir"/$(echo "$core_name" | sed 's/_bt2//')_bbdk_summary.txt \
            $bbduk_opts

    # Files that do not match the required extension
    elif [[ ! "$file" == *.@(fq|fastq|fq.gz|fastq.gz) ]]; then
        echo -e "$(basename -- "$file") will not be processed as is not a .fastq or .fq.gz file."
    fi
done

# Compress output
cd "$out_dir"
gzip *.fq

echo "Done."
echo "You can now use these 16S clean reads to:"
echo "- Perform amplicon-based analysis in Qiime2 or Mothur"
