#!/bin/bash
# Script to generate a read-based coverage table with CoverM.
# Written by: Phagomica Group
# Last updated on: 2021-02-14

##------------------------------Checking the input---------------------------##:
set -e
SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat<<HELP_USAGE
Generate read-based coverage tables with CoverM.
Usage: metabiome coverm [options] -i <in_dir> -o <out_dir> -opts coverm_options

Required:
  -i  in_dir        Directory containing paired-end reads and contigs in fastq
                    and fasta format, respectively. (The filenames of the contigs and
                    their respective paired-end reads must match)
  -o out_dir        Directory in which results will be saved. This directory
                    will be created if it doesn't exist.

Options
  -m method         Method to generate the read-based coverage table.(default=metabat).
                    (Refer to CoverM manual in order to use other methods)
  -t NUM            Number of threads. (default=1)
  -opts OPTIONS     CoverM's options.
  -h, --help        Show this help.
HELP_USAGE
}

##--------------------------Exiting if input files are missing---------------##:
validate_arguments "$#"
##----------------------Saving input orders into variables------------------##:
while (("$#")); do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )       input_dir=$(readlink -f "$2"); shift 2 ;;
        -o )       out_dir=$(readlink -m "$2"); shift 2 ;;
        -t )       threads="$2"; shift 2 ;;
        -m )       method="$2"; shift 2 ;;
        -opts )    shift; cov_opts="$@"; break ;;
        * )        echo "Option '$1' not recognized"; exit 1 ;;
    esac
done
##----------------Verify that input directory exists------------------------##:
validate_input_dir
# Create output directory if it doesn't exists
validate_output_dir
##----------------Activate conda environment--------------------------------##:
activate_env metabiome-concoct
##---------------Output info-------------------------------------------------##:
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir}"
echo "Number of threads: ${threads:=1}"
echo "Coverage method: ${method:=metabat}"
echo "Coverm version: $(coverm -V)"
echo "Coverm called with options: $cov_opts"
##--------------------------Generate the read-base coverage table file---------------------------##:
cd "$out_dir"
for file in "$input_dir"/*;do
    if [[ "$file" == *@(*_R1_*|*_1).@(fq|fastq|fq.gz|fastq.gz) ]];then
            forward_file="$file"
            core_name=$(get_core_name "$forward_file")
            contig=$(get_genome_format "$input_dir"/"$core_name")
            coverm contig -1 "$forward_file" -2 $(forward_to_reverse "$forward_file") \
                --reference "$contig" -m "$method" -t "$threads" $cov_opts > "$core_name".tsv
    fi
done
