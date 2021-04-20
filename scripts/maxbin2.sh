#!/bin/bash
# MaxBin2 wrapper script for binning contigs.
# Written by: Phagomica Group
# Last updated on:  2021-03-24

set -e
SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat<<HELP_USAGE
Generate bins from metagenomic samples with MaxBin2.
Usage: metabiome maxbin2 [options] -i <in_dir> -o <out_dir> -opts maxbin2_options

Required
  -i  in_dir        Directory containing paired-end reads and contigs in fastq and
                    fasta format, respectively.(The filenames of the contigs and
                    their respective paired-end reads must match)
  -o out_dir        Directory in which results will be saved. This directory
                    will be created if it doesn't exist.

Options:
  -a                Directory containing abundance files; each abundance filename
                    must match their respective contig filename. (if not set,
                    MaxBin2 will generate the abundance files)
  -t NUM            Number of threads. (default=1)
  -opts OPTIONS     MaxBin2's options.
  -h, --help        Show this help.
HELP_USAGE
}

# Exit if command is called with no arguments
validate_arguments "$#"

# Parse command line arguments
while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift 2 ;;
        -o )        out_dir=$(readlink -m "$2"); shift 2 ;;
        -a )        abundance_dir=$(readlink -f "$2"); shift 2 ;;
        -t )        threads="$2"; shift 2 ;;
     -opts )        shift; maxbin2_opts="$@"; break ;;
         * )        echo "Option '$1' not recognized"; exit 1 ;;
    esac
done

# Verify that input directory exists
validate_input_dir

# Create output directory if it doesn't exists
validate_output_dir

# Activate conda environment
activate_env metabiome-maxbin2

# Output info
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir}"
echo "Number of threads: ${threads:=1}"
echo "MaxBin2 version: $(run_MaxBin.pl -v)"
echo "MaxBin2 called with options: $maxbin2_opts"

# MaxBin2 binning
# if the abundance file is not provided
if [[ ! -d "$abundance_dir" ]]; then
    for file in "$input_dir"/*; do
        if [[ "$file" == *@(*_R1_*|*_1).@(fq|fastq|fq.gz|fastq.gz) ]]; then
            forward_file="$file"
            core_name=$(get_core_name "$forward_file")
            contig=$(get_genome_format "$input_dir"/"$core_name")
            create_dir "$out_dir" "$core_name" && cd "$out_dir"/"$core_name"

            # Run Maxbin2
            run_MaxBin.pl -contig "$contig" -reads "$forward_file" \
                -reads2  $(forward_to_reverse "$forward_file") -out "$core_name" \
                -thread "$threads" $maxbin2_opts
        fi
    done

# if the abundance file is provided
elif [[ -d "$abundance_dir" ]]; then
    for file in "$input_dir"/*; do
        if [[ "$file" == *.@(fna|fasta|fa) ]]; then
            contig="$file"
            core_name=$(get_core_name "$contig")
            create_dir "$out_dir" "$core_name" && cd "$out_dir"/"$core_name"

            # Run Maxbin2
            run_MaxBin.pl -contig "$contig" -abund "$abundance_dir"/*"$core_name"* \
                -out "$core_name" -thread "$threads" $maxbin2_opts
        fi
    done
fi
