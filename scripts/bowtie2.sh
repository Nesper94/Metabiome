#!/bin/bash
# Bowtie2 wrapper script for the filtering of contaminating reads from metagenomic samples:
# Written by: Phagomica Group
# Last updated on: 2020-13-11

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: metabiome bowtie2 [options] -i <in_dir> -o <out_dir>"
    echo
    echo "Required: "
    echo "  -i in_dir               Input directory containing FASTQ files."
    echo "  -o out_dir              Directory in which results will be saved. This directory"
    echo "                          will be created if it doesn't exist."
    echo
    echo "Options:"
    echo "  -ho host                Host reference genome in FASTA format."
    echo "  -ph PhiX                PhiX-174 phage reference genome in FASTA format."
    echo "  -hu human               Human reference genome in FASTA format."
    echo "  -t  NUM                 Number of threads to use. (default=1)"
    echo "  -idx index_OPTIONS      bowtie2 index builder's options."
    echo "  -opts bowtie2_OPTIONS   bowtie2's options."
    echo "  -h, --help              Show this help"
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
        * )         echo "Option '$1' not recognized"; exit 1 ;;
    esac
done

##----------------------Activate Conda environment-----------------------------##:
activate_env metabiome-preprocessing

##-------------Download Human and PhiX reference genomes-----------------##:
# Highlight: Next code downloads PhiX and Human genome and checks if the
# downloads were successfull.

for i in {1..10}; do
    if [ -e "$Human" ]; then
        echo "Human Reference Genome already downloaded"
        break
    else
        echo "Downloading Human Reference Genome"
        esearch -db nucleotide -query "GCA_000001405.28" | \
        efetch -format fasta > Human_GCA_000001405.28.fasta
        if [[ -s Human_GCA_000001405.28.fasta ]]; then
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
        efetch -format fasta > PhiX_NC_001422.1.fasta
        if [[ -s PhiX_NC_001422.1.fasta ]]; then
            echo "PhiX genome was downloaded"
            PhiX=$(readlink -f PhiX_NC_001422.1.fasta)
            break
        else
            echo "Try again downloading the PhiX genome"
        fi
    fi
done

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists
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
# First checks if the index is already generated, otherwise it will be generated.
if [ -f "$host" ]; then # checks if the host sequence file is provided.
    cat "$host" "$PhiX" "$Human" > Mixed.fasta
else
    cat "$PhiX" "$Human" > Mixed.fasta
fi
##-----------------------Index the mixed fasta-------------------------------##:
dir=$(pwd) ##current directory
[ ! -f "$dir"/Mix.3.bt2 ] \
&& { echo "Index needs to be generated:";
    bowtie2-build Mixed.fasta Mix --threads "$threads" "$index_opts";}

##-----------------------Alignment------------------------------------##:
for file in "$input_dir"/*; do
    # Paired end reads
    if [[ "$file" == @(*_R1_*|*_1).@(fq|fastq|fq.gz|fastq.gz) ]]; then
        forward_file="$file"
        core_name=$(get_core_name "$forward_file")
        bowtie2 -x Mix -1 "$forward_file" -2 $(forward_to_reverse "$forward_file") \
        --un-conc-gz "$out_dir"/$(echo "$core_name" | sed 's/_trim/_bt2/').fq.gz \
        -q -p "$threads" 2> "$out_dir"/$(echo "$core_name" | sed 's/_trim/_bt2/')_summary.txt \
        $bowtie2_opts \
        > /dev/null # Bowtie2 output to terminal is excesive and we do not need it in this case

    # Single end reads
    elif [[ ! "$file" ==  *_@(R1_*|1.|R2_*|2.)* ]] && [[ "$file" == *.@(fq|fastq|fq.gz|fastq.gz) ]]; then
        unpaired_file="$file"
        core_name=$(get_core_name "$unpaired_file")
        bowtie2 -x Mix -U "$unpaired_file" \
        --un-gz "$out_dir"/$(basename -- "$unpaired_file" | sed 's/_trim/_bt2/') \
        -q -p "$threads" 2> "$out_dir"/$(get_core_name "$unpaired_file" | sed 's/_trim/_bt2/')_summary.txt \
        $bowtie2_opts \
        > /dev/null

    # Files that do not match the required extension
    elif [[ ! "$file" == *.@(fq|fastq|fq.gz|fastq.gz) ]]; then
        echo -e "$(basename -- "$file") will not be processed as is not a .fastq or .fq.gz file."
    fi
done

# Rename bowtie2 output files
cd "$out_dir"
rename "s/.fq.1.gz/_1.fq.gz/" *.fq.1.gz 2> /dev/null
rename "s/.fq.2.gz/_2.fq.gz/" *.fq.2.gz 2> /dev/null

# Set the correct file extension (.gz) for unpaired output files
rename -f "s/fq/fq.gz/" *.fq 2> /dev/null
rename -f "s/fastq/fastq.gz/"  *.fastq 2> /dev/null

echo "Done."
echo "You can now use clean reads to:"
echo "- Assemble genomes using metaspades.sh or megahit.sh."
echo "- Perform taxonomic binning with kraken2.sh or MetaPhlAn3.sh"
echo "- Perform functional profiling using humann2.sh"
echo "- Extract 16S sequences with BBduk.sh"
