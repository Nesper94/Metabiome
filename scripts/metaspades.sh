#!/bin/bash
# Metagenomic assembly using metaSPAdes from SPAdes-3.12.0
# Written by: Phagomica_club

set -e

function usage() {
    echo "Usage: $0 -i <input directory> -o <output directory> [-t <threads>] [OPTIONS]"
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
        * )         ms_opts="$@"
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

# SPAdes #

# wget http://cab.spbu.ru/files/release3.12.0/SPAdes-3.12.0-Linux.tar.gz
# tar -xzf SPAdes-3.12.0-Linux.tar.gz #extraer archivo
# cd SPAdes-3.12.0-Linux/bin/
# export PATH="directory_path:$PATH" #add SPAdes installation directory to the PATH variable

# Activate conda environment
source activate assembly

# Run metaSPAdes #

for i in "$input_dir"/*.1.fq.gz;do

	echo "Performing PE assembly with files $(basename $i) and $(basename $i | sed 's/.1.gz/.2.gz/')"
	spades.py --meta \
    -o "$out_dir" \
    -1 "$i" `# Forward files` \
    -2 $(echo "$i" | sed 's/.1.fq.gz/.2.fq.gz/') `# Reverse sequences` \
    -t "$threads" \
    "$ms_opts:=''" # Obtain other options for metaSPAdes

    #TO DO
    #"$ms_opts:=''" # Obtain other options for metaSPAdes


done
# metaspades.sh -i input -o output -t 20
