#!/bin/bash
# Bowtie2 wrapper script for the filtering of contaminating reads from metagenomic samples:
# Written by: Phagomica Group
# Last updated on: 2020-08-19

set -e

echo "Performing bowtie2 alignment. You should have a conda environment \
in order to run this script:"

function usage() {
    echo "Usage: $0 -i <input directory> -o <output directory> \
		-ho <host_reference> -e <conda_env> [-t <threads>] -ph <PhiX_sequence \
		-hu <Human_genome>"
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
        -ho )       host=$(readlink -f "$2")
            shift
            ;;
        -e )        conda_env="$2"
            shift
            ;;
				-t )        threads="$2"
		        shift
		        ;;
				-ph )        Phix="$2"
				    shift
				    ;;
				-hu )        Human="$2"
				    shift
				    ;;

        * )         echo "Option '$1' not recognized"; exit 1
            ;;
    esac
    shift
done


# Output info
echo "Input directory: ${reads_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=4}"
echo "Conda environment: ${conda_env:?'=conda environment not set'}"
echo "Phix Genome: ${PhiX:?'=PhiX genome not set'}"
echo "Human Genome: ${Human:?'=Human genome not set'}"
echo "Kraken2 version: $(bowtie2 -v)"

if [[ ! -d "$out_dir" ]]; then  # Create output directory if it doesn't exists.
    mkdir "$out_dir"
fi

##------------------Installing required packages----------------------------##:
:'
conda activate $conda_env
env_path="$CONDA_PREFIX/bin"
echo "Here lies the packages from your environment: $env_path"
cd $env_path

if [ -e bowtie2* ];then
	echo "bowtie2 installed"
else
	echo "Installing Bowtie2" &&  conda install bowtie2 --yes
fi

if [ -e entrez-direct* ]; then
	echo "entrez-direct installed"
else
	echo "Installing entrez-direct: " && conda install entrez-direct --yes
fi'

cd $read_dir
##------------Downloading Human and PhiX reference genomes-----------------##:
:'
for i in {1..10};do
	echo "PhiX and human genomes:"
	esearch -db nucleotide -query "NC_001422.1" | \
	efetch -format fasta > NC_001422.1.fasta
	esearch -db nucleotide -query "GCA_000001405.28"  | \
	efetch -format fasta > GCA_000001405.28.fasta
	echo "Confirming downloaded genomes: "
	if  [ -e NC_001422.1.fasta ];then
			echo "PhiX genome was downloaded"
			if  [ -e GCA_000001405.28.fasta ];then
				echo "Human genome was downloaded"
				break
			fi
	fi;done
'
##-------------Concatenating genomes to be aligned-------------------------##:

#cat $host $Phage $Human  > Mixed.fasta

##----------Building genome index and bowtie alignment----------------------##:
echo "Building genome index:"
bowtie2-build Mixed.fasta Mix --threads $threads && echo 'Indexing genomes \
to filter out...'

##--------------------------Pair end (PE) alignmen--------------------------##:

	for i in $(ls *f-paired.fq.gz -1);do
		echo "PE alignment: "
		bowtie2 -x Mix -1 $i -2 $(echo $i | sed 's/f-paired/r-paired/') \
		--un-conc-gz $out_dir/$(echo $i | sed 's/f-paired.fq.gz/paired_unaligned.fq.gz/') \
		-q -p $threads --met-file $out_dir/$(echo $i | sed 's/f-paired.fq.gz/PE_summary.txt/');done

##-------------------------Single end (SE) alignment------------------------##:
	for i in $(ls *unpaired.fq.gz -1);do
		echo "SE alignment: "
		bowtie2 -x Mix -U $i \
		--un-gz $out_dir/$(echo $i | sed 's/unpaired.fastq.gz/unpaired_unaligned.fq.gz/') \
		-q -p $threads --met-file $out_dir/$(echo $i | sed 's/unpaired.fastq.gz/SE_summary.txt/');done
