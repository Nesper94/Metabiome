#!/bin/bash
# Metagenomic assembly using metaSPAdes from SPAdes-3.12.0
# Written by: Phagomica_club

set -e

function usage() {
    echo "Usage: $0 -i <input directory> -o <output directory> [-t <threads>] [MetaSPADES_OPTIONS]"
    echo "WARNING: Make sure to enclose MetaSPADES_OPTIONS within quotation marks."
    echo "Output directory will be created if it doesn't exists."
}

if [[ "$#" == 0 ]]; then
    echo "Error: No arguments given." >&2
    usage
    exit 1
fi

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
        * )         MetaSPADES_opts="$@"
            ;;
    esac
    shift
done

echo "Number of threads: ${threads:=4}"

# Verify that input directory exists
if [ ! -d "$input_dir" ]; then
   echo "$0: Error: $input_dir is not a valid directory." >&2
   exit 1
fi

# Create output directory if it doesn't exists.
if [[ ! -d "$out_dir" ]]; then
    mkdir "$out_dir"
fi

# Activate conda environment
source activate assembly

# Run metaSPAdes #

forward_file_suffix=1_paired_bt2.fq.gz
reverse_file_suffix=2_paired_bt2.fq.gz

for forward_file in "$input_dir"/*$forward_file_suffix; do

	echo "Performing PE assembly with files $(basename $forward_file) and $(basename $forward_file | sed 's/_1_bt2.fq.gz/_2_bt2.fq.gz/')"
	spades.py --meta \
    -o "$out_dir" \
    -1 "$forward_file" \
    -2 $(echo "$forward_file" | sed 's/$forward_file_suffix/$reverse_file_suffix') `# Reverse sequences` \
    -t "$threads" \
    "$MetaSPADES_opts:=''" # Obtain other options for metaSPAdes

done
