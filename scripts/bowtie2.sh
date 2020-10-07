#!/bin/bash
# Bowtie2 wrapper script for the filtering of contaminating reads from metagenomic samples:
# Written by: Phagomica Group
# Last updated on: 2020-10-06

set -e

function usage() {
    echo "Usage: $0 -i <in dir> -o <out dir> -ho <host> -ph <PhiX> -hu <human> [-e <conda_env>] [-t <threads>]"
    echo ""
    echo "<in dir>  Input directory containing FASTQ files."
    echo "<out dir> Directory in which results will be saved. This directory"
    echo "          will be created if it doesn't exists."
    echo "<host>    Host reference genome in FASTA format."
    echo "<PhiX>    PhiX-174 phage reference genome in FASTA format."
    echo "<human>   Human reference genome in FASTA format."
    echo ""
    echo "Options:"
    echo "-t        Number of threads to use."
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
        -ho )       host=$(readlink -f "$2")
            shift
            ;;
        -e )        conda_env="$2"
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

# Verify that input directory exists
if [ ! -d "$input_dir" ]; then
   echo "$0: Error: $input_dir is not a valid directory."
   exit 1
fi

if [[ ! -d "$out_dir" ]]; then  # Create output directory if it doesn't exists.
    mkdir "$out_dir"
fi

# Activate conda environment
source activate preprocessing

# Output info
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=4}"
echo "Conda environment: ${conda_env:?'=conda environment not set'}"
echo "Phix Genome: ${PhiX:?'=PhiX genome not set'}"
echo "Human Genome: ${Human:?'=Human genome not set'}"
echo "Bowtie2 version: $(bowtie2  --version)"

##--------------Download Human and PhiX reference genomes-----------------##:

for i in {1..10};do
	echo "PhiX and human genomes:"
	esearch -db nucleotide -query "NC_001422.1" |	efetch -format fasta > NC_001422.1.fasta \
	esearch -db nucleotide -query "GCA_000001405.28" | efetch -format fasta > GCA_000001405.28.fasta \
	echo "Confirming downloaded genomes: "
	if  [ -e NC_001422.1.fasta ];then
			echo "PhiX genome was downloaded"
			if  [ -e GCA_000001405.28.fasta ];then
				echo "Human genome was downloaded"
				break
			fi
	fi;done

##-----------------Concatenate genomes to be aligned------------------------##:

cat "$host" "$PhiX" "$Human"  > Mixed.fasta

##---------------Build genome index and bowtie alignment-----------------##:
echo "Building genome index:"
bowtie2-build Mixed.fasta Mix --threads $threads && echo 'Indexing genomes \
to filter out...'

##-----------------Perform Paired end reads (PE) alignment-------------------##:

for i in "$input_dir"/*1_paired_trim.fq.gz;do
    
    echo "Paired End alignment: "
    
    bowtie2 -x Mix -1 "$i" \
    -2 $(echo $i | sed 's/1_paired_trim/2_paired_trim/') \
    --un-conc-gz "$out_dir"/$(echo $(basename -- $i) | sed 's/1_paired_trim.fq.gz/paired_bt2.fq.gz/') \
    -q -p $threads 2> "$out_dir"/$(echo $(basename -- $i) | sed 's/1_paired_trim.fq.gz/paired_bt2_summary.txt/')
    done
    

##-----------------Perform Single end reads (SE) alignment-------------------##:

for i in "$input_dir"/*unpaired_trim.fq.gz;do
    
    echo "Single End alignment: "
    
    bowtie2 -x Mix -U "$i" \
    --un-gz "$out_dir"/$(echo $(basename -- $i) | sed 's/unpaired_trim.fq.gz/unpaired_bt2.fq.gz/') \
    -q -p $threads 2> "$out_dir"/$(echo $(basename -- $i) | sed 's/unpaired_trim.fq.gz/unpaired_bt2_summary.txt/')
    done

echo "Done."
echo "You can now use clean reads to:"
echo "- Assemble genomes using metaspades.sh or megahit.sh."
echo "- Perform taxonomic binning with kraken2.sh or MetaPhlAn3.sh"
echo "- Perform functional profiling using humann2.sh"
echo "- Extract 16S sequences with BBduk.sh"
