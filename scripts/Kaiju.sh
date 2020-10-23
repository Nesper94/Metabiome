#!/bin/bash
# Kaiju wrapper script for the taxonomic binning from metagenomic samples
# Written by: Phagomica Group
# Last updated on: 2020-08-19

set -e

function usage() {
    echo "Usage: $0 -i <input directory> -o <output directory> \
		-D <database_directory> -d <name_kaiju_database> -e <conda_env> [-t <threads>]"
    echo ""
    echo "<input_directory>  Input directory containing FASTQ files."
    echo "<output_directory> Directory in which results will be saved. This directory"
    echo "          will be created if it doesn't exists."
    echo "<database_directory>    Kaiju Database directory."
    echo "<name_kaiju_database>   Kaiju Database name. Please refer to the Kaiju"
    echo "github repository in order to know the names"
    echo ""
    echo "Options:"
    echo "<threads>        Number of threads to use."
}

if [[ "$#" == 0 ]]; then
    echo "No arguments given." >&2
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
        -D )        database=$(readlink -f "$2")
            shift
            ;;
        -d )        name_kaijudb="$2"
            shift
            ;;
        -t )        threads="$2"
           shift
            ;;
        * )        echo "Option '$1' not recognized"; exit 1
            ;;
    esac
    shift
done

# Verify that input directory exists
if [ ! -d "$input_dir" ]; then
  echo "$0: Error: $input_dir is not a valid directory." >&2
  exit 1
fi
# Create output directory if it doesn't exists.
if [[ ! -d "$out_dir" ]]; then
    mkdir "$out_dir"
fi

# Output info
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir:?'Output directory not set'}"
echo "Number of threads: ${threads:=4}"
echo "Kaiju reference database: ${database:?'=Kaiju database directory not set'}"
echo "Kaiju database name: ${name_kaijudb:?'=Kaiju database name not set'}"
echo "Kaiju version: $(kaiju -h)"

# Activate conda environment
source activate binning

###----------------------Making Kaiju Database------------------------------###:
cd "$database"
kaiju-makedb -s "$name_kaijudb" -t "$threads"

###----------------------Paired-end reads-----------------------------------###:
for i in "$input_dir"/*paired_unaligned.fq.1.gz;do
  kaiju -t "$database"/*nodes.dmp -f "$database"/"$name_kaijudb"/*.fmi -i "$i" \
  -j $(echo $i | sed 's/.1.gz/.2.gz/') -o  "$out_dir"/$(echo $(basename -- $i) | sed 's/.fq.1.gz/_kaiju_out.txt/') \
  -z "$threads";done

###----------------------Single-end reads-----------------------------------###:
for i in "$input_dir"/*unpaired_unaligned.fq.gz;do
  kaiju -t "$database"/*nodes.dmp -f "$database"/"$name_kaijudb"/*.fmi -i "$i" \
  -o "$out_dir"/$(echo $(basename -- $i) | sed 's/.fq.gz/_kaiju_out.txt/') -z "$threads";done


###----------------------Convention applied-------------------------------###:

  1_bt2.fq.gz
  2_bt2.fq.gz
  un_bt2.fq.gz

  ###----------------------Paired-end reads-----------------------------------###:
  for i in "$input_dir"/*1_bt2.fq.gz;do
    kaiju -t "$database"/*nodes.dmp -f "$database"/"$name_kaijudb"/*.fmi -i "$i" \
    -j $(echo $i | sed 's/_1_/_2_/') \
    -o  "$out_dir"/$(echo $(basename -- $i) | sed 's/1_bt2.fq.gz/pe_kaiju.txt/') \
    -z "$threads";done

  ###----------------------Single-end reads-----------------------------------###:
  for i in "$input_dir"/*un_bt2.fq.gz;do
    kaiju -t "$database"/*nodes.dmp -f "$database"/"$name_kaijudb"/*.fmi -i "$i" \
    -o "$out_dir"/$(echo $(basename -- $i) | sed 's/bt2.fq.gz/se_kaiju.txt/') -z "$threads";done



###----------------------Kaiju to Krona-------------------------------------###:
echo "conda install Krona"

if [ -d "$krona" ]; then
    for i in "$out_dir"/*kaiju_out.txt;do
          kaiju2krona -t "$database"/nodes.dmp  -n "$database"/names.dmp \
          -i $i -o "$krona"/$(echo $(basename -- $i) | sed 's/kaiju_out.txt/krona.txt/');done
    for i in "$krona"/*.txt;do
      ktImportText $i -o "$html"/$(echo $(basename -- $i) | sed 's/.txt/.html/');done
fi
###----------------------Classification summary from output-------------------------------------###:
if [ -d "$class" ]; then
  for i in "$out_dir"/*kaiju_out.txt;do
        kaiju2table -t "$database"/nodes.dmp  -n "$database"/names.dmp \
        -p -m 1.0 -o "$class"/$(echo $(basename -- $i) | sed 's/kaiju_out.txt/classif.txt/');done
fi

###---------------------Adding taxa names to output file----------------------------------------###:
if [ -d "$taxa_names" ]; then
      for i in "$out_dir"/*kaiju_out.txt;do
                    kaiju-addTaxonNames -t "$database"/nodes.dmp -p \
                    -n "$database"/names.dmp -i $i \
                    -o "$taxa_names"/$(echo $(basename -- $i) | sed 's/kaiju_out.txt/tax_names.txt/');done
fi

###-----------------------------Merging outputs from kraken and Kaiju--------------------------###:
if [ -d "$merge" ] && [ -d "$class" ]; then
    for i in "$class"/*classif.txt;do
          for j in "$kraken"/*.txt;do
                  kaiju-mergeOutputs  -i <(sort -k2,2 $i) -j <(sort -k2,2 $j) \
                  -t "$database"/nodes.dmp -c lowest \
                  -o "$taxa_names"/$(echo $(basename -- $i) | sed 's/kaiju_out.txt/classif.txt/') -v ;done
     ;done
fi
