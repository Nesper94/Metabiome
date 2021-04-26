#!/bin/bash
# Kaiju wrapper script for the taxonomic binning from metagenomic samples
# Written by: Phagomica Group
# Last updated on: 2021-01-25

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat <<HELP_USAGE
Generate taxonomic bins of metagenomic samples with kaiju.
Usage: metabiome kaiju [options] -i <in_dir> -o <out_dir> -D <db_dir> -d <kaiju_db> -opts kaiju_options

Required:
  -i in_dir         Input directory containing FASTQ reads files.
  -o out_dir        Directory in which results will be saved. This directory
                    will be created if it doesn't exist.
  -D db_dir         Directory where to download Kaiju Database.
  -d kaiju_db       Kaiju Database name. Please refer to the Kaiju github repository.

Options:
  -t NUM            Number of threads to use (default=1).
  -c class_dir      Output directory of the kaiju taxa classification.
  -x taxa_dir       Output directory of the kaiju taxa name assignation.
  -m merge_dir      Output directory of the kaiju and kraken2 merged files.
  -k krona_dir      Output directory of the kaiju-to-krona figures.
  -kr kraken2_dir   Directory containing Kraken2 output files.(Must be set if merge option is set.)
  -opts OPTIONS     Kaiju's options.
  -h, --help        Show this help.
HELP_USAGE
}
##-----------------------Exit if called with no arguments--------------------##:
validate_arguments "$#"
##------------------Saving input orders into variables-----------------------##:
while (("$#")); do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift 2 ;;
        -o )        out_dir=$(readlink -m "$2"); shift 2 ;;
        -D )        database=$(readlink -f "$2"); shift 2 ;;
        -d )        name_kaijudb="$2"; shift 2 ;;
        -t )        threads="$2"; shift 2 ;;
        -c )        class=$(readlink -f "$2"); shift 2 ;;
        -x )        taxa_names=$(readlink -f "$2"); shift 2 ;;
        -m )        merge=$(readlink -f "$2"); shift 2 ;;
        -k )        krona=$(readlink -f "$2"); shift 2 ;;
       -kr )        kraken2=$(readlink -f "$2"); shift 2 ;;
     -opts )        shift; kaiju_opts="$@"; break ;;
        * )         echo "Option '$1' not recognized"; exit 1 ;;
    esac
done

##----------------Verify that input directory exists-------------------------##:
validate_input_dir
##---------------Create output directory if it doesn't exists----------------##:
validate_output_dir
##------------------------Activate conda environment-------------------------##:
activate_env metabiome-taxonomic-binning
##------------------------------Output info----------------------------------##:
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: $input_dir"
echo "Output directory: $out_dir"
echo "Number of threads: ${threads:=1}"
echo "Kaiju reference database: ${database:?'=Kaiju database directory not set'}"
echo "Kaiju database name: ${name_kaijudb:?'=Kaiju database name not set'}"
echo "Kaiju called with options: $kaiju_opts"
##---------------------------Make Kaiju Database-----------------------------##:
if [[ ! -d "$database"/"$name_kaijudb" ]];then
    echo "Database $name_kaijudb needs to be generated:"
    cd "$database" && kaiju-makedb -s "$name_kaijudb" -t "$threads"
fi
##-----------------------Taxonomic binning------------------------------------##:
  for file in "$input_dir"/*; do
    # Paired end reads
    if [[ "$file" == @(*_R1_*|*_1).@(fq|fastq|fq.gz|fastq.gz) ]]; then
        forward_file="$file"
        core_name=$(get_core_name "$forward_file")
        kaiju -t "$database"/*nodes.dmp -f "$database"/"$name_kaijudb"/*.fmi -i "$forward_file" \
            -j $(forward_to_reverse "$forward_file") \
            -o  "$out_dir"/$(get_core_name "$forward_file" | sed 's/_bt2/_kaiju/').txt \
            -z "$threads" $kaiju_opts

    # Single end reads
    elif [[ ! "$file" ==  *_@(*R1_*|*1.|*R2_*|*2.)* ]] && [[ "$file" == *.@(fq|fastq|fq.gz|fastq.gz) ]]; then
        unpaired_file="$file"
        core_name=$(get_core_name "$unpaired_file")
        kaiju -t "$database"/*nodes.dmp -f "$database"/"$name_kaijudb"/*.fmi -i "$unpaired_file" \
            -o "$out_dir"/$(get_core_name "$unpaired_file" | sed 's/_bt2/_kaiju/').txt \
            -z "$threads" $kaiju_opts

    #Files that do not match the required extension:
    elif [[ ! "$file" == *.@(fq|fastq|fq.gz|fastq.gz) ]]; then
        echo -e "$(basename -- "$file") will not be processed as is not a .fastq or .fq.gz file."
    fi
 done
##-----------------------------Kaiju to Krona--------------------------------##:
if [ -d "$krona" ]; then
    cd "$krona" && mkdir html && cd ..
    for kaiju_out in "$out_dir"/*.txt;do
        kaiju2krona -t "$database"/nodes.dmp  -n "$database"/names.dmp \
            -i "$kaiju_out" -o "$krona"/$(get_core_name "$kaiju_out" | sed 's/_kaiju/.krona/')
    done
    for krona_out in "$krona"/*krona*;do
        ktImportText "$krona_out" -o "$krona"/html/$(get_core_name "$krona_out" | sed 's/.krona/_krona/').html
    done
fi
##---------------Classification summary from Kaiju's output------------------##:
if [ -d "$class" ]; then
    for kaiju_out in "$out_dir"/*.txt;do
        kaiju2table "$kaiju_out" -t "$database"/nodes.dmp  -n "$database"/names.dmp \
            -p -m 1.0 -r genus -o "$class"/$(get_core_name "$kaiju_out" | sed 's/kaiju/classif/').txt
    done
fi
##-----------------------Add taxa names to output file-----------------------##:
if [ -d "$taxa_names" ]; then
    for kaiju_out in "$out_dir"/*.txt;do
        kaiju-addTaxonNames -t "$database"/nodes.dmp -p \
            -n "$database"/names.dmp -i "$kaiju_out" \
            -o "$taxa_names"/$(get_core_name "$kaiju_out" | sed 's/kaiju/tax_names/').txt
    done
fi
##-------------------Merge outputs from Kraken and Kaiju---------------------##:
if [ -d "$merge" ] && [ -d "$kraken2" ]; then
    cd "$merge"
    for kaiju_out in "$out_dir"/*_paired_*;do
        core_name=$(get_core_name "$kaiju_out")
        echo "sorting  kaiju output files"
        sort -k2,2 "$kaiju_out" > $(echo "$core_name" | sed 's/kaiju/kaiju.mrgsrt/')
        sort -k2,2 "$kraken2"/$(echo "$core_name" | sed 's/kaiju/kraken2_out.tsv/') > \
            $(echo "$core_name" | sed 's/kaiju/kraken2.mrgsrt/')
    done
    for sorted_file in "$merge"/*kaiju.mrgsrt*;do
        kaiju-mergeOutputs  -i "$sorted_file" \
            -j $(basename -- "$sorted_file" | sed 's/kaiju/kraken2/') \
            -t "$database"/nodes.dmp -c lowest \
            -o "$merge"/$(get_core_name "$sorted_file" | sed 's/kaiju/merged/').txt -v
    done
    rm *mrgsrt*
echo "You can re-run the classification, taxonomic and krona analysis on these merged files,"
fi
