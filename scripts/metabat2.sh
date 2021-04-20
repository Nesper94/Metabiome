#!/bin/bash
# MetaBAT2 wrapper script for the binning of contig assemblies.
# Written by: Phagomica Group
# Last updated on: 2021-03-24

set -e
SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat<<HELP_USAGE
Generate bins from metagenomic samples with MetaBAT2.
Usage: metabiome metabat2 [options] -i <in_dir> -o <out_dir> -opts metabat2_options

Required:
  -i  in_dir        Directory containing contigs in gzipped fasta format.
  -o out_dir        Directory in which results will be saved. This directory
                    will be created if it doesn't exists.

Options:
  -co cov_dir       Directory containing MetaBAT2's reads coverage file; each
                    coverage filename must match to their respective contig
                    filename.(These coverage files can be generated with
                    'metabiome coverm')
  -t NUM            Number of threads. (default=1)
  -opts OPTIONS     MetaBAT2's options.
  -h, --help        Show this help.
HELP_USAGE
}

# Exit if command is called with no arguments
validate_arguments "$#"

# Parse command line arguments
while (("$#")); do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )       input_dir=$(readlink -f "$2"); shift 2 ;;
        -o )       out_dir=$(readlink -m "$2"); shift 2 ;;
        -co )      cov_dir=$(readlink -f "$2"); shift 2 ;;
        -t )       threads="$2"; shift 2 ;;
        -opts )    shift; metabat2_opts="$@"; break ;;
        * )        echo "Option '$1' not recognized"; exit 1 ;;
    esac
done

# Verify that input directory exists
validate_input_dir

# Create output directory if it doesn't exists
validate_output_dir

# Activate conda environment
activate_env metabiome-metabat2

# Output info
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir}"
echo "Number of threads: ${threads:=1}"
echo "MetaBAT2 version: $(metabat2 -h)"
echo "MetaBAT2 called with options: $metabat2_opts"

# MetaBAT2 binning
# if the read-based coverage file is not provided
if [[ ! -d "$cov_dir" ]]; then
    for file in "$input_dir"/*; do
        if [[ "$file" == *.@(fna.gz|fasta.gz|fa.gz) ]]; then
            contig="$file"
            core_name=$(get_core_name "$contig")
            create_dir "$out_dir" "$core_name"

            # run MetaBAT2
            metabat2 -i "$contig" \
                -o "$out_dir"/"$core_name"/"$core_name" \
                -t "$threads" $metabat2_opts
        fi
    done

# if the read-based coverage file is provided
elif [[ -d "$cov_dir" ]]; then
    for file in "$input_dir"/*; do
        if [[ "$file" == *.@(fna.gz|fasta.gz|fa.gz) ]]; then
            contig="$file"
            core_name=$(get_core_name "$contig")
            create_dir "$out_dir" "$core_name"

            # run MetaBAT2
            metabat2 -i "$contig" \
                --abdFile "$cov_dir"/*"$core_name"* \
                -o "$out_dir"/"$core_name"/"$core_name" \
                -t "$threads" $metabat2_opts
        fi
    done
fi
