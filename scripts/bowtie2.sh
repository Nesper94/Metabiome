#!/bin/bash
# Bowtie2 wrapper script for the filtering of contaminating reads from metagenomic samples:
# Written by: Phagomica Group
# Last updated on: 2020-13-11

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: metabiome bowtie2 [options] -i <in_dir> -o <out_dir>"
    echo ""
    echo "Mandatory: "
    echo "-i, <in_dir>               Input directory containing FASTQ files."
    echo "-o, <out_dir>              Directory in which results will be saved. This directory"
    echo "                           will be created if it doesn't exist."
    echo "Options:"
    echo "-ho, <host>                Host reference genome in FASTA format. (optional)"
    echo "-ph, <PhiX>                PhiX-174 phage reference genome in FASTA format. (optional)"
    echo "-hu, <human>               Human reference genome in FASTA format. (optional)"
    echo "-t,  <threads>             Number of threads to use. (default=1)"
    echo "-idx, <index_OPTIONS>      bowtie2 index builder's options (optional)."
    echo "-opts, <bowtie2_OPTIONS>   bowtie2's options (optional)."
    echo "Please execute '(bowtie2 --help)' in order to see bowtie2's documentation "

}

# Exit if command is called with no arguments
validate_arguments "$#"

##---------------------Save input parameters into variables------------------##:

while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift 2 ;;
        -o )        out_dir=$(readlink -m "$2"); shift 2 ;;
        -ho )       host=$(readlink -f "$2"); shift 2 ;;
        -t )        threads="$2"; shift 2 ;;
        -ph )       PhiX="$2"; shift 2 ;;
        -hu )       Human="$2"; shift 2 ;;
        -idx )      shift; index_opts="$@";;
        -opts )     shift; bowtie2_opts="$@"; break ;;
        * )        echo "Option '$1' not recognized"; exit 1 ;;
    esac
done

##----------------------Activate Conda environment-----------------------------##:
activate_env metabiome-preprocessing

##------------Download Human and PhiX reference genomes-----------------##:
##Highlight: Next code downloads PhiX and Human genome and checks if the downloads
##were successfull.

for i in {1..10}; do
    if [ -e "$Human" ]; then
      echo "Human Reference Genome already downloaded"
      break
    else
      echo "Downloading Human Reference Genome"
      esearch -db nucleotide -query "NC_001422.1" | \
      efetch -format fasta > Human_GCA_000001405.28.fasta
      if  [[ -e Human_GCA_000001405.28.fasta ]] && [[ -s Human_GCA_000001405.28.fasta ]]; then
          echo "Human genome was downloaded"
          Human=$(readlink -f Human_GCA_000001405.28.fasta)
          break
      else
        echo "Try again downloading the human genome"
      fi
    fi
done

for i in {1..10}; do
  if [ -e "$PhiX" ]; then
    echo "PhiX Genome already downloaded"
    break
  else
    echo "Downloading PhiX Reference Genome"
    esearch -db nucleotide -query "NC_001422.1" | \
    efetch -format fasta > NC_001422.1.fasta
    if  [[ -e NC_001422.1.fasta ]] && [[ -s NC_001422.1.fasta ]]; then
      echo "PhiX genome was downloaded"
      PhiX=$(readlink -f NC_001422.1.fasta)
      break
    else
      echo "Try again downloading the PhiX genome"
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
echo "Phix Genome: ${PhiX:?'=PhiX genome not set'}"
echo "Human Genome: ${Human:?'=Human genome not set'}"
echo "Bowtie2 called with options: $bowtie2_opts"


##------------Concatenate genomes to be aligned and build genome index-------##:
##First checks if the index is already generated, otherwise it will be generated.
if [ -f "$host" ];then ##checks if the host sequence file is provided.
    cat "$host" "$PhiX" "$Human" > Mixed.fasta
else
    cat "$PhiX" "$Human" > Mixed.fasta
fi
##--------------------Indexing the mixed fasta-------------------------------##:
dir=$(pwd) ##current directory
[ ! -f "$dir"/Mix.3.bt2 ] \
&& { echo "Index needs to be generated:";
    bowtie2-build Mixed.fasta Mix --threads "$threads" "$index_opts";}

##--------------------------Pair end (PE) alignment--------------------------##:

for forward_file in "$input_dir"/*1_paired_* ; do
    echo "Performing paired reads alignment..."
    bowtie2 -x Mix -1 "$forward_file" -2 $(echo "$forward_file" | sed 's/_1_paired_/_2_paired_/') \
    --un-conc-gz "$out_dir"/$(echo $(basename -- "$forward_file" ) | sed 's/1_paired_trim/paired_bt2/') \
    -q -p "$threads" 2> "$out_dir"/$(echo $(basename -- "$forward_file" ) | sed 's/1_paired_trim.fq.gz/paired_bt2_summary.txt/') \
    "$bowtie2_opts" \
    > /dev/null # Bowtie2 output to terminal is excesive and we do not need it in this case
done

##-------------------------Single end (SE) alignment------------------------##:
for unpaired_file in "$input_dir"/*_unpaired_* ; do
    echo "Performing single reads alignment..."
    bowtie2 -x Mix -U "$unpaired_file" --un-gz "$out_dir"/$(echo $(basename -- "$unpaired_file" ) | sed 's/unpaired_trim/unpaired_bt2/') \
    -q -p "$threads" 2> "$out_dir"/$(echo $(basename -- "$unpaired_file" ) | sed 's/unpaired_trim.fq.gz/unpaired_bt2_summary.txt/') \
    "$bowtie2_opts" \
    > /dev/null
done

##---------------- Rename files according to the naming convention-----------##:
cd "$out_dir"
##renaming paired-end reads:
for i in {1..2}; do
  rename 's/paired_bt2.fq.'"$i"'.gz/'"$i"'_paired_bt2.fq.gz/' *paired_bt2.fq."$i".gz
done


echo "Done."
echo "You can now use clean reads to:"
echo "- Assemble genomes using metaspades.sh or megahit.sh."
echo "- Perform taxonomic binning with kraken2.sh or MetaPhlAn3.sh"
echo "- Perform functional profiling using humann2.sh"
echo "- Extract 16S sequences with BBduk.sh"
