#!/bin/bash
# Metagenomic assembly using MEGAHIT-1.2.9
# Written by: Phagomica_club

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: metabiome megahit -i <input directory> -o <output directory> [-t <threads>] [OPTIONS]"
    echo "  Output directory will be created if it doesn't exists."
    echo ""
    echo "This is a wrapper script that passes OPTIONS to MEGAHIT, so this is the"
    echo "documentation for MEGAHIT:"
    echo ""
    echo "$(megahit --help)"
}

# Exit if command is called with no arguments
validate_arguments "$#"

# Activate conda environment
activate_env metabiome-genome-assembly

while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0
            ;;
        -i )        input_dir=$(readlink -f "$2")
            shift
            ;;
        -o )        out_dir=$(readlink -f "$2")
            shift
            ;;
        -t )        threads="$2"
            shift
            ;;
        --k-list)   kmer_list="$2"
		shift
            ;;
        * )         MEGAHIT_opts="$@"
            ;;
    esac
    shift
done

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists.
validate_output_dir

#Output info
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=1}"
echo "MEGAHIT version: $(megahit --version)"

# Run MEGAHIT on paired-end files

forward_file_suffix=1_paired_bt2.fq.gz
reverse_file_suffix=2_paired_bt2.fq.gz

for forward_file in "$input_dir"/*$forward_file_suffix; do

	echo "Performing PE assembly with $(basename $forward_file) and $(basename "$forward_file" | sed "s/$forward_file_suffix/$reverse_file_suffix/")"
    out_name=$(basename "$forward_file" | sed 's/$forward_file_suffix/.fq.gz/')
	megahit \
    -o "$out_dir"/"$out_name" \
    -1 "$forward_file" `#Forward files(1 files, paired with files in "$pe2")` \
    -2 $(echo "$forward_file" | sed "s/$forward_file_suffix/$reverse_file_suffix/") `# Reverse files (2 files, paired with files in "$pe1"` \
    -t "$threads" "${MEGAHIT_opts:=}" \
    --presets meta-large # Optimization for large & complex metagenomes, like soil

done
