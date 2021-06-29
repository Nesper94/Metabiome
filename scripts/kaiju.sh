#!/bin/bash
# Kaiju wrapper script for the taxonomic binning of metagenomic samples.
# Written by: Phagomica Group
# Last updated on: 2021-05-24

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat <<HELP_USAGE
Generate taxonomic bins of metagenomic samples with Kaiju.
Usage: metabiome kaiju [options] -i <in_dir> -o <out_dir> -d <kaiju_db> -opts kaiju_options

Required:
  -i in_dir         Input directory containing FASTQ files.
  -o out_dir        Directory in which results will be saved. This directory will be created if it
                    doesn't exist.

  -d kaiju_db       Kaiju Database name. (Please refer to the Kaiju GitHub repository)

Options:
  -t NUM            Number of threads to use (default=1).
  -D db_dir         Directory containing Kaiju's Database. If this directory does not exist, it will
                    be created and the database downloaded automatically.
  -c                Create classification summary of Kaiju's output files.
  -x                Add taxa names to Kaiju's output files.
  -k                Generate krona graphs of Kaiju's output files.
  -opts OPTIONS     Kaiju's options.
  -h, --help        Show this help.
  -hh               Show Kaiju's help message.
HELP_USAGE
}

# Exit if called with no arguments
validate_arguments "$#"

# Get input parameters
while (("$#")); do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -hh )       activate_env metabiome-taxonomic-binning; kaiju -h; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift 2 ;;
        -o )        out_dir=$(readlink -m "$2"); shift 2 ;;
        -D )        database=$(readlink -m "$2"); shift 2 ;;
        -d )        name_kaijudb="$2"; shift 2 ;;
        -t )        threads="$2"; shift 2 ;;
        -c )        class=class; shift 1 ;;
        -x )        taxa_names=taxa_names; shift 1 ;;
        -k )        krona=krona; shift 1 ;;
     -opts )        shift; kaiju_opts="$@"; break ;;
        * )         echo "Option '$1' not recognized"; exit 1 ;;
    esac
done

# Verify that input directory exists
validate_input_dir

# Create output directory if it doesn't exists
validate_output_dir

# Activate conda environment
activate_env metabiome-taxonomic-binning

# Output info
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: $input_dir"
echo "Output directory: $out_dir"
echo "Number of threads: ${threads:=1}"
echo "Kaiju database name: ${name_kaijudb:?'Kaiju database name not set'}"
echo "Kaiju called with options: $kaiju_opts"

# Make Kaiju Database
# Check if the database directory is provided
if [[ ! -d "$database" ]]; then
        create_dir "$out_dir" kaiju_db
        database="$out_dir"/kaiju_db
        echo "the database will be created in $database"
        cd "$database"

        # Make kaiju database
        kaiju-makedb -s "$name_kaijudb" -t "$threads"

# If the directory is provided but the specific database is not downloaded
elif [[ -d "$database" ]] && [[ ! -d "$database"/"$name_kaijudb" ]]; then
        echo "Database $name_kaijudb needs to be generated:"
        cd "$database"

        # Make kaiju database
        kaiju-makedb -s "$name_kaijudb" -t "$threads"
fi

# Taxonomic binning
for file in "$input_dir"/*; do
    # Paired end reads
    if [[ "$file" == @(*_R1_*|*_1).@(fq|fastq|fq.gz|fastq.gz) ]]; then
        forward_file="$file"
        core_name=$(get_core_name "$forward_file")
        kaiju -t "$database"/*nodes.dmp \
            -f "$database"/"$name_kaijudb"/*.fmi \
            -i "$forward_file" \
            -j $(forward_to_reverse "$forward_file") \
            -o  "$out_dir"/$(get_core_name "$forward_file" | sed 's/_bt2/_kaiju/').txt \
            -z "$threads" $kaiju_opts

    # Unpaired reads
    elif [[ ! "$file" ==  *_@(*R1_*|*1.|*R2_*|*2.)* ]] && [[ "$file" == *.@(fq|fastq|fq.gz|fastq.gz) ]]; then
        unpaired_file="$file"
        core_name=$(get_core_name "$unpaired_file")
        kaiju -t "$database"/*nodes.dmp \
            -f "$database"/"$name_kaijudb"/*.fmi \
            -i "$unpaired_file" \
            -o "$out_dir"/$(get_core_name "$unpaired_file" | sed 's/_bt2/_kaiju/').txt \
            -z "$threads" $kaiju_opts

    # Files that do not match the required extension:
    elif [[ ! "$file" == *.@(fq|fastq|fq.gz|fastq.gz) ]]; then
        echo -e "$(basename -- "$file") will not be processed as is not a .fastq or .fq.gz file."
    fi
done

# Kaiju to Krona
if [[ -n "$krona" ]]; then
    # Create output directory for krona figures
    create_dir "$out_dir" "$krona" && create_dir "$out_dir"/"$krona" html
    cd "$out_dir"/"$krona"

    # Generate krona figures from kaiju output files
    for kaiju_out in "$out_dir"/*.txt; do
        kaiju2krona -t "$database"/nodes.dmp \
            -n "$database"/names.dmp \
            -i "$kaiju_out" \
            -o "$out_dir"/"$krona"/$(get_core_name "$kaiju_out" | sed 's/_kaiju/.krona/')
    done

    # Generate the .html of the krona figures
    for krona_out in "$out_dir"/"$krona"/*krona*; do
        ktImportText "$krona_out" \
            -o "$out_dir"/"$krona"/html/$(get_core_name "$krona_out" | sed 's/.krona/_krona/').html
    done
fi

# Classification summary of Kaiju's output
if [[ -n "$class" ]]; then
    # Create output directory for classification summary
    create_dir "$out_dir" "$class"

    # Generate classification summary of Kaiju output files
    for kaiju_out in "$out_dir"/*.txt; do
        kaiju2table "$kaiju_out" \
            -t "$database"/nodes.dmp  \
            -n "$database"/names.dmp \
            -p -m 1.0 -r genus \
            -o "$out_dir"/"$class"/$(get_core_name "$kaiju_out" | sed 's/kaiju/classif/').txt
    done
fi

# Add taxa names to output file
if [[ -n "$taxa_names" ]]; then
    # Create output directory for the assignation of taxa names
    create_dir "$out_dir" "$taxa_names"

    # Generate assignation of taxa names of Kaiju output files
    for kaiju_out in "$out_dir"/*.txt; do
        kaiju-addTaxonNames -t "$database"/nodes.dmp -p \
            -n "$database"/names.dmp \
            -i "$kaiju_out" \
            -o "$out_dir"/"$taxa_names"/$(get_core_name "$kaiju_out" | sed 's/kaiju/tax_names/').txt
    done
fi
