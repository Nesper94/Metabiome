#!/bin/bash
# BBduk wrapper script for the DNAr 16S picking strategy from metagenomic samples:
# Written by: Phagomica Group
# Last updated on: 2020-08-19

set -e

echo "Performing kmer 16S rDNA kmer picking strategy with BBDuk: "

function usage() {
    echo "Usage: $0 -i <input directory> -o <output directory> \
		-D <16S_DATABASE> [-e <conda_env>] [-t <threads>]"
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
        -i )        input_dir=$(readlink -f "$2")
            shift
            ;;
        -o )        out_dir=$(readlink -f "$2")
            shift
            ;;
        -D )        database=$(readlink -f "$2")
            shift
            ;;
        -e )        conda_env="$2"
            shift
            ;;
        -t )        threads="$2"
           shift
            ;;
        * )        echo "Option '$1' not recognized"; exit 1
            ;;
    esac
    shift
done


# Output info
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=4}"
echo "Conda environment: ${conda_env:?'=conda environment not set'}"
echo "Reference database: ${database:?'=reference database not set'}"

# Verify that input directory exists
if [ ! -d "$input_dir" ]; then
   echo "$0: Error: $input_dir is not a valid directory."
   exit 1
fi

if [[ ! -d "$out_dir" ]]; then  # Create output directory if it doesn't exists.
    mkdir "$out_dir"
fi

###---------------Moving to your conda environment packages location---------##:
echo 'Lets activate your environment: ' && conda activate "$conda_env"
##Path to your conda environment
echo "Here lies the packages from your environment: "
env_path= echo $CONDA_PREFIX/bin
cd "$env_path"

###############Installing required packages###################:

if [ -e bbmap* ];then
	echo "bbmap is installed"
else
	echo "Installing bbmap" &&  conda install bbmap --yes
fi

cd "$input_dir"

##Matching reads to the 16S rDNA SSU from SILVA Database,

#For PE reads:



for i in "$input_dir"/*paired_unaligned.fq.1.gz;do
  bbduk.sh in=$i in2=$(echo $i | 's/.1.gz/.2.gz/') \
	ref="$database" outm="$out_dir"/$(echo $(basename -- $i) | sed 's/.1.gz/BBduk_1_16S.fq/') \
  outm2="$out_dir"/$(echo $(basename -- $i) | sed 's/.1.gz/BBduk_2_16S.fq/') \
	outs="$out_dir"/$(echo $(basename -- $i) | sed 's/.1.gz/BBduk_single_16S.fq/') \
  stats="$out_dir"/$(echo $(basename -- $i) | sed 's/fq.1.gz/BBduk_stats.txt/') ordered=T ;done

#For SE reads:
for i in "$input_dir"/*unpaired_unaligned.fq.gz;do
  bbduk.sh in=$i ref="$database" outm="$out_dir"/$(echo $(basename -- $i) | sed 's/.fq.gz/-16S.fq/') \
	stats="$out_dir"/$(echo $(basename -- $i) | sed 's/unpaired_unaligned.fq.gz/unpaired_statistics_16S.txt/') \
  ordered=T;done

####Compressing################:
cd $out_dir
gzip *.fq
