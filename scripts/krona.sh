#!/bin/bash
# Krona wrapper script
# Authored by: Estefany Lorenzana López, Cristian Grisales Vargas and
# Juan Camilo Arboleda Rivera

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: metabiome krona [Options] -i <input directory> -o <output directory>"
    echo
    echo "Required:"
    echo "  -i in_dir   Input directory containing Kraken2 output files."
    echo "  -o out_dir  Directory in which results will be saved. This directory"
    echo "              will be created if it doesn't exist."
    echo "Options:"
    echo "  -h, --help  Show this help"
}

# Exit if command is called with no arguments
validate_arguments "$#"

while (("$#")); do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift 2 ;;
        -o )        out_dir=$(readlink -m "$2"); shift 2 ;;
        * )         echo "Unknown option: $1"; exit 1 ;;
    esac
done

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists.
validate_output_dir

# Activate conda environment
activate_env metabiome-taxonomic-binning

# Output info
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: $input_dir"
echo "Output directory: $out_dir"
echo

if [[ ! -f "$CONDA_PREFIX"/opt/krona/taxonomy/taxonomy.tab ]]; then
    echo "Updating Krona taxonomy..."
    ktUpdateTaxonomy.sh "$CONDA_PREFIX"/opt/krona/taxonomy/
fi

cd "$out_dir"

for kraken2_out in "$input_dir"/*_kraken2_out.tsv; do
    krona_input=$(basename -- "$kraken2_out" | sed "s/_kraken2_out\.tsv/\.krona/")
    cat "$kraken2_out" | cut -f 2,3 > "$krona_input"
done

ktImportTaxonomy *.krona

rm *.krona # Clean up intermediary files
