#!/bin/bash
# Metagenomic assembly using MEGAHIT-1.2.9
# Written by: Phagomica_club

set -e

function usage() {
    echo "Usage:"
    echo "  $0 -i <input directory> -o <output directory> [-t <threads>] [OPTIONS]"
    echo "  Output directory will be created if it doesn't exists."
    echo ""
    echo "This is a wrapper script that passes OPTIONS to MEGAHIT, so this is the"
    echo "documentation for MEGAHIT:"
    echo ""
    echo "$(megahit --help)"
}

if [[ "$#" == 0 ]]; then
    echo "No arguments given." >&2
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
        --k-list)   kmer_list="$2"
		shift
            ;;
        * )         MEGAHIT_opts="$@"
            ;;
    esac
    shift
done

# Activate conda environment
source activate assembly

#Output info
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=4}"
echo "MEGAHIT version: $(megahit --version)"

# Verify that input directory exists
if [ ! -d "$input_dir" ]; then
   echo "$0: Error: $input_dir is not a valid directory." >&2
   exit 1
fi

# Create output directory if it doesn't exists.
if [[ ! -d "$out_dir" ]]; then
    mkdir "$out_dir"
fi

# Run MEGAHIT on paired-end files

forward_file_suffix=1_paired_bt2.fq.gz
reverse_file_suffix=2_paired_bt2.fq.gz

for forward_file in "$input_dir"/*$forward_file_suffix; do

	echo "Performing PE assembly with $(basename $forward_file) and $(basename "$forward_file" | sed 's/f-paired/r-paired/')"
    
    out_name=$(basename "$forward_file" | sed 's/1_paired_bt2//')
    
	megahit \
    -o "$out_dir"/"$out_name" \
    -1 "$forward_file" `#Forward files(1 files, paired with files in "$pe2")` \
    -2 $(echo "$forward_file" | sed 's/$forward_file_suffix/$reverse_file_suffix/') `# Reverse files (2 files, paired with files in "$pe1"` \
    -t "$threads" "$MEGAHIT_opts:=''" \
    --presets meta-large # Optimization for large & complex metagenomes, like soil

done
