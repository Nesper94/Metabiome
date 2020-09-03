#!/bin/bash
# Metagenomic assembly using MEGAHIT-1.2.9
# Written by: Phagomica_club

sudo apt update -y
set -e


function usage() {
    echo "Usage: $0 -i <input directory> -o <output directory> [-t <threads>] [-k <kmers>]"
    echo "Output directory will be created if it doesn't exists."
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
        -i )        reads_dir=$(readlink -f "$2")
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

# Download and install #
wget https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
tar zvxf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
cd MEGAHIT-1.2.9-Linux-x86_64-static/bin/
#./megahit --test  # run on a toy dataset
#./megahit -1 MY_PE_READ_1.fq.gz -2 MY_PE_READ_2.fq.gz -o MY_OUTPUT_DIR

#MEGAHIT options

for i in "$reads_dir"/*f-paired*.fq.gz;do
	echo "PE assembly"
	megahit
    -o "$out_dir" \
    -1 "$i" `#Forward files(1 files, paired with files in "$pe2")` \
    -2 $(echo "$i" | sed 's/.1.gz/.2.gz/') `# Reverse files (2 files, paired with files in "$pe1"` \
    -t/--num-cpu-threads "$2"
    --no-mercy
    --k-list "21,33,55" `# list of k-mer sizes to be use separeted comma(no more than 141)`

done

#--k-min "$" minimum kmer size (<= {0}), must be odd number [21]
#--k-max "$#"
