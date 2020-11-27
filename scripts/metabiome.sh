#!/bin/bash
# ---------------------
# Main Metabiome script
# ---------------------

set -e

VERSION=1.0.0
SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh
echo "SCRIPTS_DIR"
function usage() {
    echo "Usage: metabiome [Commands|Options]"
    echo ""
    echo "  Commands:"
    echo "    fastqc        Check read quality with FastQC"
    echo "    trimmomatic"
    echo "    bowtie2       Remove contaminant sequences"
    echo "    metaspades    Assemble reads into contigs"
    echo "    megahit       Assemble reads into contigs"
    echo ""
    echo "  Options:"
    echo "    -h, --help  Show this help"
}

# Exit if command is called with no arguments
validate_arguments "$#"

# Get input arguments
if [[ "$1" = "bbduk" ]]; then
  "$SCRIPTS_DIR"/BBduk.sh ${@:2}

elif [[ "$1" = "bowtie2" ]]; then
  "$SCRIPTS_DIR"/bowtie2.sh ${@:2}

elif [[ "$1" = "humann2" ]]; then
  "$SCRIPTS_DIR"/humann2.sh ${@:2}

elif [[ "$1" = "kaiju" ]]; then
  "$SCRIPTS_DIR"/Kaiju.sh ${@:2}

elif [[ "$1" = "kraken2" ]]; then
  "$SCRIPTS_DIR"/kraken2.sh ${@:2}

elif [[ "$1" = "megahit" ]]; then
  "$SCRIPTS_DIR"/megahit.sh ${@:2}

elif [[ "$1" = "metaphlan3" ]]; then
  "$SCRIPTS_DIR"/MetaPhlAn3.sh ${@:2}

elif [[ "$1" = "metaspades" ]]; then
  "$SCRIPTS_DIR"/metaspades.sh ${@:2}

elif [[ "$1" = "trimmomatic" ]]; then
  "$SCRIPTS_DIR"/trimmomatic.sh ${@:2}

elif [[ "$1" = "-h" ]] || [[ "$1" = "--help" ]]; then
  usage

elif [[ "$1" = "-v" ]] || [[ "$1" = "--version" ]]; then
  echo "Metabiome $VERSION"
fi
