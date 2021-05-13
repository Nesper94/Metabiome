#!/bin/bash
# Wrapper script to refine bins with DAS tool.
# Written by: Phagomica Group
# Last updated on: 2021-02-14

set -e
SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat <<HELP_USAGE
Calculate non-redundant and optimized bins with DAS Tool.
Usage: metabiome das_tool [options] -i <in_dir> -o <out_dir> -opts DAS_Tool_options

Required:
-i in_dir           Input directory containing scaffolds-to-bins tsv files and
                    contigs in fasta format. (Scaffolds-to-bins tsv filenames must
                    match to their respective contigs filenames)
-o out_dir          Directory in which results will be saved. This directory
                    will be created if it doesn't exist.

Options:
-pt proteins_dir    Directory containing the predicted proteins in prodigal fasta
                    format (.faa) of each contig. (The predicted protein fasta filenames
                    must match to their respective contig filenames)
-t  NUM             Number of threads to use. (default=1)
-opts OPTIONS       DAS Tool's options.
-h, --help          Show this help.

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
        -t )       threads="$2"; shift 2 ;;
        -pt )      proteins=$(readlink -m "$2"); shift 2 ;;
        -opts )    shift; das_opts="$@"; break ;;
        * )        echo "Option '$1' not recognized"; exit 1 ;;
    esac
done

# Verify that input directory exists
validate_input_dir

# Create output directory if it doesn't exists
validate_output_dir

# Activate conda environment
activate_env metabiome-das_tool

# Output info
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Number of threads: ${threads:=1}"
echo "$(DAS_Tool -v)"
echo "DAS tool called with options: $das_opts"

# Calculate bins with DAS Tool

# if the predicted protein files are not provided
if [[ ! -d "$proteins" ]]; then
    for file  in "$input_dir"/*; do
        if [[ "$file" == *.@(fa|fna|fasta) ]]; then
            contig="$file"
            core_name=$(get_core_name "$contig")
            create_dir "$out_dir" "$core_name"


            # Generate a comma-separated list of the scaffolds-to-bins tsv files
            tsv_list=$(ls -1 "$input_dir"/*"$core_name"*.tsv | paste -d, -s)

            # Run DAS_Tool
            DAS_Tool --bins "$tsv_list" \
                -c "$contig" \
                -o "$out_dir"/"$core_name"/"$core_name" \
                --threads "$threads" \
                $das_opts


        elif [[ ! "$file" == *.@(fa|fna|fasta) ]] && [[ ! "$file" == *.@(tsv) ]]; then
            echo "$file is not either a bin or a scaffold-to-bin tsv file, so it will be ignore"
        fi
    done

# if the predicted protein files are provided
elif [[ -d "$proteins" ]]; then
    for file  in "$input_dir"/*; do
        if [[ "$file" == *.@(fa|fna|fasta) ]]; then
            contig="$file"
            core_name=$(get_core_name "$contig")
            create_dir "$out_dir" "$core_name"

            # Generate a comma-separated list of the scaffolds-to-bins tsv files
            tsv_list=$(ls -1 "$input_dir"/*"$core_name"*.tsv | paste -d, -s)


            # Prodigal fasta file
            prodigal_file=$(echo "$proteins"/*"$core_name"*.faa)

            # Check if the prodigal fasta file is correctly provided
            if [[ -e "$prodigal_file" ]]; then

                # Run DAS Tool
                echo "Skipping Prodigal gene prediction"
                DAS_Tool --bins "$tsv_list" \
                    -c "$contig" \
                    -o "$out_dir"/"$core_name"/"$core_name" \
                    --threads "$threads" \
                    --proteins "$prodigal_file" \
                    $das_opts


            elif [[ ! -e "$prodigal_file" ]]; then
                echo "There is not a Prodigal gene prediction file ($core_name*.faa) in $proteins"
                exit 1
            fi


        elif [[ ! "$file" == *.@(fa|fna|fasta) ]] && [[ ! "$file" == *.@(tsv) ]]; then
            echo "$file is not either a bin or a scaffold-to-bin tsv file, so it will be ignore"
        fi
    done
fi
