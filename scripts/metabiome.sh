#!/bin/bash
# ---------------------
# Main Metabiome script
# ---------------------

set -e

VERSION=1.0.0
SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat<<HELP_USAGE
Usage: metabiome [Commands|Options]

  Preprocessing:
    bowtie2       Remove contaminant sequences.
    qc            Check read quality with FastQC and MultiQC.
    trimmomatic   Perform quality trimming on Illumina sequence data.

  Taxonomic binning:
    kaiju         Perform taxonomic classification of sequences.
    kraken2       Perform taxonomic classification of sequences.
    krona         Create Krona charts using Kraken2 output.

  Taxonomic profiling:
    humann        Profile the abundance of microbial metabolic pathways and molecular functions.
    metaphlan3    Generate taxonomic profiling bins.

  Extract 16S rDNA sequences:
    bbduk         Extract 16S rDNA sequences.

  Genome Assembly:
    megahit       Assemble reads into contigs.
    metaquast     Evaluate genome assembly.
    metaspades    Assemble reads into contigs.

  Genome binning:
    concoct       Generate bins from contigs and paired-end reads.
    coverm        Create read-coverage tables.
    maxbin2       Generate bins from contigs and paired-end reads.
    metabat2      Generate bins from contigs and paired-end reads.

  Bins assessment:
    das_tool      Calculate non-redundant and optimized bins.

  Options:
    -h, --help    Show this help.
    -v, --version Show Metabiome's version.

HELP_USAGE
}

# Exit if command is called with no arguments
validate_arguments "$#"

# Get input arguments
if [[ "$1" = "bbduk" ]]; then
  bash "$SCRIPTS_DIR"/bbduk.sh ${@:2}

elif [[ "$1" = "bowtie2" ]]; then
  bash "$SCRIPTS_DIR"/bowtie2.sh ${@:2}

elif [[ "$1" = "humann" ]]; then
  bash "$SCRIPTS_DIR"/humann3.sh ${@:2}

elif [[ "$1" = "kaiju" ]]; then
  bash "$SCRIPTS_DIR"/kaiju.sh ${@:2}

elif [[ "$1" = "kraken2" ]]; then
  bash "$SCRIPTS_DIR"/kraken2.sh ${@:2}

elif [[ "$1" = "krona" ]]; then
  "$SCRIPTS_DIR"/krona.sh ${@:2}

elif [[ "$1" = "megahit" ]]; then
  bash "$SCRIPTS_DIR"/megahit.sh ${@:2}

elif [[ "$1" = "metaphlan3" ]]; then
  bash "$SCRIPTS_DIR"/metaphlan3.sh ${@:2}

elif [[ "$1" = "metaquast" ]]; then
  bash "$SCRIPTS_DIR"/metaquast.sh ${@:2}

elif [[ "$1" = "metaspades" ]]; then
  bash "$SCRIPTS_DIR"/metaspades.sh ${@:2}

elif [[ "$1" = "qc" ]]; then
  bash "$SCRIPTS_DIR"/qc.sh ${@:2}

elif [[ "$1" = "trimmomatic" ]]; then
  bash "$SCRIPTS_DIR"/trimmomatic.sh ${@:2}

elif [[ "$1" = "metabat2" ]]; then
  bash "$SCRIPTS_DIR"/metabat2.sh ${@:2}

elif [[ "$1" = "maxbin2" ]]; then
  bash "$SCRIPTS_DIR"/maxbin2.sh ${@:2}

elif [[ "$1" = "concoct" ]]; then
  bash "$SCRIPTS_DIR"/concoct.sh ${@:2}

elif [[ "$1" = "coverm" ]]; then
  bash "$SCRIPTS_DIR"/coverm.sh ${@:2}

elif [[ "$1" = "das_tool" ]]; then
  bash "$SCRIPTS_DIR"/das_tool.sh ${@:2}

elif [[ "$1" = "-h" ]] || [[ "$1" = "--help" ]]; then
  usage

elif [[ "$1" = "-v" ]] || [[ "$1" = "--version" ]]; then
  echo "Metabiome $VERSION"

else
  echo "Error: Option $1 not recognized."
  exit 1
fi
