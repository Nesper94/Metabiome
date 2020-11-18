#!/bin/bash
# Kaiju wrapper script for the taxonomic binning from metagenomic samples
# Written by: Phagomica Group
# Last updated on: 2020-10-28

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: metabiome kaiju -i <input directory> -o <output directory> -D <database_directory> -d <name_kaiju_database> [-t <threads>] [-k krona_dir] [-c class] [-x taxa_names] [-m merge ]"
    echo ""
    echo "Options:"
    echo "<input_directory>  Input directory containing clean FASTQ files."
    echo "<output_directory> Directory in which results will be saved. This directory"
    echo "will be created if it doesn't exist."
    echo "<database_directory>    Kaiju Database directory."
    echo "<name_kaiju_database>   Kaiju Database name. Please refer to the Kaiju"
    echo "github repository in order to know the names"
    echo "<threads> Number of threads to use."
    echo "<class>        Directory where to put the classification summary of Kaiju's output. (optional)"
    echo "<taxa_names>   Directory where to put the taxa names from the Kaiju's output file. (optional)"
    echo "<merge>        Directory where to put the merging outputs from Kaiju and Kraken. Needs classification (class) output in order to run (optional)"
    echo "<krona_dir>    Directory where to put Krona output. (optional)"
    echo "<kraken_input> Kraken output files's directory. (optional)"
}

##--------------------------Exiting if input files are missing---------------##:
validate_arguments "$#"

##------------------Saving input orders into variables-----------------------##:
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
        -D )        database=$(readlink -f "$2")
            shift
            ;;
        -d )        name_kaijudb=$(readlink -f "$2")
            shift
            ;;
        -t )        class=$(readlink -f "$2")
           shift
            ;;
        -t )        taxa_names=$(readlink -f "$2")
           shift
            ;;
        -t )        merge=$(readlink -f "$2")
           shift
            ;;
        -t )        krona=$(readlink -f "$2")
           shift
            ;;
        -k )        kraken=$(readlink -f "$2")
           shift
            ;;
        * )        echo "Option '$1' not recognized"; exit 1
            ;;
    esac
    shift
done

##----------------Verify that input directory exists-------------------------##:
validate_input_dir

##---------------Create output directory if it doesn't exists----------------##:
validate_output_dir


##---------------------Activate conda environment----------------------------##:
activate_env binning

##---------------Output info-------------------------------------------------##:
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=1}"
echo "Kaiju reference database: ${database:?'=Kaiju database directory not set'}"
echo "Kaiju database name: ${name_kaijudb:?'=Kaiju database name not set'}"
echo "Kaiju version: $(kaiju -h)"

###----------------------Making Kaiju Database-------------------------------##:
cd "$database"
kaiju-makedb -s "$name_kaijudb" -t "$threads"

##-----------------------Paired-end reads------------------------------------##:
for i in "$input_dir"/*1_paired_bt2.fq.gz;do
  kaiju -t "$database"/*nodes.dmp -f "$database"/"$name_kaijudb"/*.fmi -i "$i" \
  -j $(echo $i | sed 's/1_paired/2_paired/') \
  -o  "$out_dir"/$(echo $(basename -- $i) | sed 's/1_paired_bt2.fq.gz/paired_kaiju.txt/') \
  -z "$threads";done

##-----------------------Single-end reads------------------------------------##:
for i in "$input_dir"/*unpaired_bt2.fq.gz;do
  kaiju -t "$database"/*nodes.dmp -f "$database"/"$name_kaijudb"/*.fmi -i "$i" \
  -o "$out_dir"/$(echo $(basename -- $i) | sed 's/unpaired_bt2.fq.gz/unpaired_kaiju.txt/') \
  -z "$threads";done

##-----------------------Kaiju to Krona--------------------------------------##:

if [ -d "$krona" ]; then
    mkdir html
    for i in "$out_dir"/*kaiju.txt;do
          kaiju2krona -t "$database"/nodes.dmp  -n "$database"/names.dmp \
          -i $i -o "$krona"/$(echo $(basename -- $i) | sed 's/kaiju.txt/krona.txt/');done
    for i in "$krona"/*.txt;do
      ktImportText $i -o "$krona"/html/$(echo $(basename -- $i) | sed 's/.txt/.html/');done
fi
##-------------------Classification summary from Kaiju's output----------------------##:
if [ -d "$class" ]; then
  for i in "$out_dir"/*kaiju.txt;do
        kaiju2table -t "$database"/nodes.dmp  -n "$database"/names.dmp \
        -p -m 1.0 -o "$class"/$(echo $(basename -- $i) | sed 's/kaiju.txt/classif.txt/');done
fi

##------------------Adding taxa names to output file-------------------------##:
if [ -d "$taxa_names" ]; then
      for i in "$out_dir"/*kaiju.txt;do
        kaiju-addTaxonNames -t "$database"/nodes.dmp -p \
        -n "$database"/names.dmp -i $i \
        -o "$taxa_names"/$(echo $(basename -- $i) | sed 's/kaiju.txt/tax_names.txt/');done
fi

##---------------Merging outputs from kraken and Kaiju-----------------------##:
if [ -d "$merge" ] && [ -d "$class" ] && [ -d "$kraken"]; then
    for i in "$class"/*classif.txt;do
          for j in "$kraken"/*.txt;do
                  kaiju-mergeOutputs  -i <(sort -k2,2 $i) -j <(sort -k2,2 $j) \
                  -t "$database"/nodes.dmp -c lowest \
                  -o "$taxa_names"/$(echo $(basename -- $i) | sed 's/kaiju.txt/classif.txt/') -v ;done
    done
fi
