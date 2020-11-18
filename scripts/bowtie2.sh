ls#!/bin/bash
# Bowtie2 wrapper script for the filtering of contaminating reads from metagenomic samples:
# Written by: Phagomica Group
# Last updated on: 2020-13-11

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: metabiome bowtie2 -i <in dir> -o <out dir> -ho <host> [-ph <PhiX>] [-hu <human>] [-t <threads>] "
    echo ""
    echo "Options:"
    echo "<in dir>  Input directory containing FASTQ files."
    echo "<out dir> Directory in which results will be saved. This directory"
    echo "will be created if it doesn't exist."
    echo "<host>    Host reference genome in FASTA format. (optional)"
    echo "<threads> Number of threads to use. (optional)"
    echo "<PhiX>    PhiX-174 phage reference genome in FASTA format. (optional)"
    echo "<human>   Human reference genome in FASTA format. (optional)"
}

# Exit if command is called with no arguments
validate_arguments "$#"

##---------------------Save input parameters into variables------------------##:
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
        -ho )       host=$(readlink -f "$2")
            shift
            ;;
        -t )        threads="$2"
           shift
            ;;
        -ph )       PhiX="$2"
           shift
            ;;
        -hu )       Human="$2"
           shift
            ;;
        * )        echo "Option '$1' not recognized"; exit 1
            ;;
    esac
    shift
done

##----------------------Activate Conda environment-----------------------------##:
activate_env preprocessing

##------------Download Human and PhiX reference genomes-----------------##:
##Highlight: Next code downloads PhiX and Human genome and checks if the downloads
##were successfull.

for i in {1..10}; do

    if [ -e "$Human" ]; then
      echo "Human Reference Genome already downloaded"
    else
      echo "Downloading Human Reference Genome"
      esearch -db nucleotide -query "GCA_000001405.28" | \
      efetch -format fasta > GCA_000001405.28.fasta
      if  [ -e GCA_000001405.28.fasta ]; then
          echo "Human genome was downloaded"
          Human=$(readlink -f GCA_000001405.28.fasta)
      fi
    fi

    if [ -e "$PhiX" ]; then
      echo "PhiX Genome already downloaded"
      break
    else
      echo "Downloading PhiX Reference Genome"
      esearch -db nucleotide -query "NC_001422.1" | \
      efetch -format fasta > NC_001422.1.fasta
      if  [ -e NC_001422.1.fasta ]; then
        echo "PhiX genome was downloaded"
        PhiX=$(readlink -f NC_001422.1.fasta)
        break
      fi
    fi
done

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists.
validate_output_dir

# Output info
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: ${input_dir}"
echo "Output directory: ${out_dir}"
echo "Number of threads: ${threads:=1}"
echo "Host Genome: ${host:?'=Host genome not set'}"
echo "Phix Genome: ${PhiX:?'=PhiX genome not set'}"
echo "Human Genome: ${Human:?'=Human genome not set'}"
echo "Bowtie2 version: $(bowtie2  --version)"

##------------Concatenate genomes to be aligned and build genome index-------##:
##First checks if the index is already generated, otherwise it will be generated.
if [ -f "$host" ];then ##checks if the host sequence file is provided.
    cat "$host" "$PhiX" "$Human" > Mixed.fasta
else
    cat "$PhiX" "$Human" > Mixed.fasta
fi
##--------------------Indexing the mixed fasta-------------------------------##:
dir=$(pwd) ##current directory
[ ! -f "$dir"/Mix.3.bt2 ] \ ##Checks if the index is already generated
&& { echo "Index needs to be generated:";
    bowtie2-build Mixed.fasta Mix --threads "$threads";}

##--------------------------Pair end (PE) alignment--------------------------##:
for i in "$input_dir"/*1_paired_trim.fq.gz; do
    echo "Performing paired reads alignment..."
    bowtie2 -x Mix -1 $i -2 $(echo $i | sed 's/1_paired_trim/2_paired_trim/') \
    --un-conc-gz "$out_dir"/$(echo $(basename -- $i) | sed 's/1_paired_trim/paired_bt2/') \
    -q -p $threads 2> "$out_dir"/$(echo $(basename -- $i) | sed 's/1_paired_trim.fq.gz/paired_bt2_summary.txt/') \
    > /dev/null # Bowtie2 output to terminal is excesive and we don't need it in this case
done

##-------------------------Single end (SE) alignment------------------------##:

for i in "$input_dir"/*unpaired_trim.fq.gz; do
    echo "Performing single reads alignment..."
    bowtie2 -x Mix -U $i --un-gz "$out_dir"/$(echo $(basename -- $i) | sed 's/unpaired_trim/unpaired_bt2/') \
    -q -p $threads 2> "$out_dir"/$(echo $(basename -- $i) | sed 's/unpaired_trim.fq.gz/unpaired_bt2_summary.txt/') \
    > /dev/null
done

##---------------- Rename files according to the naming convention-----------##:
cd "$out_dir"
rename 's/paired_bt2.fq.1.gz/1_paired_bt2.fq.gz/' *paired_bt2.fq.1.gz
rename 's/paired_bt2.fq.2.gz/2_paired_bt2.fq.gz/' *paired_bt2.fq.2.gz

echo "Done."
echo "You can now use clean reads to:"
echo "- Assemble genomes using metaspades.sh or megahit.sh."
echo "- Perform taxonomic binning with kraken2.sh or MetaPhlAn3.sh"
echo "- Perform functional profiling using humann2.sh"
echo "- Extract 16S sequences with BBduk.sh"
