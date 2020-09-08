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
    echo "No arguments given."
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
        * )         echo "Option '$1' not recognized"; exit 1
            ;;
    esac
    shift

done
#Output info
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=4}"
echo "MEGAHIT version: $(megahit --version)"

# Verify that input directory exists
if [ ! -d "$input_dir" ]; then
   echo "$0: Error: $input_dir is not a valid directory."
   exit 1
fi

# Create output directory if it doesn't exists.
if [[ ! -d "$out_dir" ]]; then
    mkdir "$out_dir"
fi

# For the moment this script only runs MEGAHIT on paired-end files

for i in "$input_dir"/*f-paired*.fq.gz;do
	echo "Performing PE assembly with $(basename $i) and $(basename "$i" | sed 's/f-paired/r-paired/')"
    out_name=$(basename "$i" | sed 's/f-paired//')
	megahit \
    -o "$out_dir"/"$out_name" \
    -1 "$i" `#Forward files(1 files, paired with files in "$pe2")` \
    -2 $(echo "$i" | sed 's/f-paired/r-paired/') `# Reverse files (2 files, paired with files in "$pe1"` \
    -t "$threads" \
    --presets meta-large # Optimization for large & complex metagenomes, like soil
#   --no-mercy `# Qué hace esto?`\
#   --k-list "21,33,55" `# list of k-mer sizes to be use separeted comma(no more than 141)`
#   -r/--read  <se>   comma-separated list of fasta/q single-end files
done
