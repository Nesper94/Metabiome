#!/bin/bash
# Kaiju wrapper script for the taxonomic binning from metagenomic samples
# Written by: Phagomica Group
# Last updated on: 2020-10-28

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: metabiome kaiju [options] -i <in_dir> -o <out_dir> -D <db_dir> -d <kaiju_db>"
    echo ""
    echo "Mandatory:"
    echo "  -i in_dir       Input directory containing FASTQ files."
    echo "  -o out_dir      Directory in which results will be saved. This directory"
    echo "                  will be created if it doesn't exist."
    echo "  -D db_dir       Kaiju Database directory."
    echo "  -d kaiju_db     Kaiju Database name. Please refer to the Kaiju github repository"
    echo
    echo "Options:"
    echo "  -t NUM          Number of threads to use (default=1)."
    echo "  -c class_dir    Directory of the classification summary of Kaiju's output."
    echo "  -x taxa_dir     Directory of the taxa names from the Kaiju's output file."
    echo "  -m merge_dir    Directory of the merging outputs from Kaiju and Kraken."
    echo "  -k krona_dir    Directory of the Krona output."
    echo "  -kr kraken_in   Kraken output files's directory."
    echo "  -opts OPTIONS   Kaiju's options."
    echo
    echo "Please execute '(kaiju -h)' in order to see Kaiju's documentation"
}

##--------------------------Exiting if input files are missing---------------##:
validate_arguments "$#"

##------------------Saving input orders into variables-----------------------##:
while [[ -n "$1" ]]; do
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
       -kr )        kraken=$(readlink -f "$2"); shift 2 ;;
     -opts )        shift; kaiju_opts="$@"; break ;;
        * )         echo "Option '$1' not recognized"; exit 1 ;;
    esac
done

##----------------Verify that input directory exists-------------------------##:
validate_input_dir

##---------------Create output directory if it doesn't exists----------------##:
validate_output_dir

##---------------------Activate conda environment----------------------------##:
activate_env metabiome-taxonomic-binning

##---------------Output info-------------------------------------------------##:
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=1}"
echo "Kaiju reference database: ${database:?'=Kaiju database directory not set'}"
echo "Kaiju database name: ${name_kaijudb:?'=Kaiju database name not set'}"
echo "Kaiju version: $(kaiju -h)"

###----------------------Making Kaiju Database-------------------------------##:
[ ! -f "$database"/"$name_kaijudb"/assembly_summary.txt ] \
&& { echo "Database $name_kaijudb needs to be generated:";
cd "$database";kaiju-makedb -s "$name_kaijudb" -t "$threads"; }

##-----------------------Paired-end reads------------------------------------##:
for forward_file in "$input_dir"/*1_paired*;do
  kaiju -t "$database"/*nodes.dmp -f "$database"/"$name_kaijudb"/*.fmi -i "$forward_file" \
  -j $(echo "$forward_file" | sed 's/1_paired/2_paired/') \
  -o  "$out_dir"/$(echo $(basename -- "$forward_file") | sed 's/1_paired_bt2.fq.gz/paired_kaiju.txt/') \
  -z "$threads" "$kaiju_opts"
done

##-----------------------Single-end reads------------------------------------##:
for forward_file in "$input_dir"/*_unpaired_*;do
  kaiju -t "$database"/*nodes.dmp -f "$database"/"$name_kaijudb"/*.fmi -i "$forward_file" \
  -o "$out_dir"/$(echo $(basename -- "$forward_file" ) | sed 's/unpaired_bt2.fq.gz/unpaired_kaiju.txt/') \
  -z "$threads" "$kaiju_opts"
done



##-----------------------Kaiju to Krona--------------------------------------##:
if [ -d "$krona" ]; then
  cd "$krona" && mkdir html && cd ..
    for kaiju_out in "$out_dir"/*kaiju.txt;do
          kaiju2krona -t "$database"/nodes.dmp  -n "$database"/names.dmp \
          -i "$kaiju_out" -o "$krona"/$(echo $(basename -- "$kaiju_out") | sed 's/_kaiju.txt/.krona/')
    done
    for krona_out in "$krona"/*.txt;do
          ktImportText "$krona_out" -o "$krona"/html/$(echo $(basename -- "$krona_out") | sed 's/.krona/krona.html/')
    done
fi
##-------------------Classification summary from Kaiju's output----------------------##:
if [ -d "$class" ]; then
  for kaiju_out in "$out_dir"/*kaiju.txt;do
        kaiju2table "$kaiju_out" -t "$database"/nodes.dmp  -n "$database"/names.dmp \
        -p -m 1.0 -r genus -o "$class"/$(echo $(basename -- "$kaiju_out") | sed 's/kaiju.txt/classif.txt/')
  done
fi
##------------------Adding taxa names to output file-------------------------##:
if [ -d "$taxa_names" ]; then
  for kaiju_out in "$out_dir"/*kaiju.txt;do
        kaiju-addTaxonNames -t "$database"/nodes.dmp -p \
        -n "$database"/names.dmp -i "$kaiju_out" \
        -o "$taxa_names"/$(echo $(basename -- "$kaiju_out") | sed 's/kaiju.txt/tax_names.txt/')
  done
fi
##---------------Merging outputs from kraken and Kaiju-----------------------##:
if [[ -d "$merge" ]] && [[ -d "$kraken"]] ; then
    for kaiju_out in "$out_dir"/*kaiju.txt;do
          for kraken_out in "$kraken"/*.fq;do
                  echo "sorting kaiju and kraken output files"
                  sort -k2,2 "$kaiju_out" > $(echo $(basename -- "$kaiju_out" ) | sed 's/kaiju.txt/kaiju.mrgsrt/')
                  sort -k2,2 "$kraken_out" > $(echo $(basename -- "$kraken_out" ) | sed 's/.fq/kraken.mrgsrt/')
                  kaiju-mergeOutputs  -i *kaiju.mrgsrt -j *kraken.mrgsrt \
                  -t "$database"/nodes.dmp -c lowest \
                  -o "$merge"/$(echo $(basename -- "$kaiju_out") | sed 's/kaiju.txt/merged.txt/') -v
                  rm *.mrgsrt
          done
    done
fi
