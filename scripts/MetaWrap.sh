#!/bin/bash
# MetaWrap binning module's wrapper script for the binning of draft assemblies.
# Written by: Phagomica Group
# Last updated on: 2020-10-28

##------------------------------Checking the input---------------------------##:
set -e
SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: metabiome metawrap -i <input_directory> -a <assembly> [-r <RAM>] [-t <threads>] [-met <metabat2_dir>] [-max <maxbin2_dir>] [-cnc <concoct_dir>] [-check <checkM>]"
    echo ""
    echo "Options:"
    echo "<input_directory>  Input directory containing clean FASTQ files."
    echo "will be created if it doesn't exist."
    echo "<assembly>  Draft metagenome assemblies."
    echo "<RAM>  Memory RAM. (optional)"
    echo "<threads>  Number of threads. (optional)"
    echo "<metabat2>  metabat2 output directory. Set only if metabat2 is one of the chosen binner to use in metawrap binning module. (optional)"
    echo "<maxbin2>  maxbin2 output directory. Set only if maxbin2 is one of the chosen binner to use in metawrap binning module. (optional)"
    echo "<concoct>  concoct2 output directory. Set only if concoct is one of the chosen binner to use in metawrap binning module. (optional)"
    echo "<checkM>  checkM output directory. Set only if checkM is chosen as the quality assessment tool of the metagenomic bins. (optional)"
}

##--------------------------Exiting if input files are missing---------------##:
validate_arguments "$#"
##----------------------Saving input orders into variables------------------##:
while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0
            ;;
        -i )        input_dir=$(readlink -f "$2")
            shift
            ;;
        -a )        assembly=$(readlink -f "$2")
            shift
            ;;
        -r )        ram="$2"
            shift
            ;;
        -t )        threads="$2"
           shift
            ;;
        -met )      metabat2=$(readlink -f "$2")
           shift
            ;;
        -max )      maxbin2=$(readlink -f "$2")
           shift
            ;;
        -cnc )      concoct=$(readlink -f "$2")
           shift
            ;;
        -check )    checkM=$(readlink -f "$2")
           shift
            ;;
        * )        echo "Option '$1' not recognized"; exit 1
            ;;
    esac
    shift
done

##----------------Verify that input directory exists------------------------##:
validate_input_dir
##---------------Create output directory if it doesn't exists---------------##:
validate_output_dir
##----------------Activate conda environment--------------------------------##:
activate_env metawrap
##---------------Output info-------------------------------------------------##:
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Assemblies directory: ${assembly:?'Assemblies directory not set'}"
echo "Number of threads: ${threads:=1}"
echo "Kaiju reference database: ${database:?'=Kaiju database directory not set'}"
echo "Kaiju database name: ${name_kaijudb:?'=Kaiju database name not set'}"
echo "Kaiju version: $(metawrap -v)"

##--------------Software to choose in metawrap binning module---------------##:

soft=()
if [ -d "$maxbin2" ];then
  echo "maxbin2 selected"
  soft+=("--maxbin2")
fi
if [ -d "$metabat2" ];then
  echo "metabat2 selected"
  soft+=("--metabat2")
fi
if [ -d "$concoct" ];then
  echo "concoct selected"
  soft+=("--concoct")
fi
if [ -d "$checkM" ];then
  echo "checkM selected"
  soft+=("--checkM")
fi

if [ ${#soft[@]} -eq 0 ];then
  echo "none of the binners from metawrap module were chosen"; exit 1
else
  echo "at least one binner was succesfully chosen"
fi

##-----------------suffixes from clean fastq input files---------------------##:
forward_file_suffix=1_paired_bt2.fq.gz
reverse_file_suffix=2_paired_bt2.fq.gz
forward_unpaired_file_suffix=1_unpaired_bt2.fq.gz
reverse_unpaired_file_suffix=2_unpaired_bt2.fq.gz

##-----------------------binning---------------------------------------------##:

for binner in "${soft[@]}";do
  for forward_file in "$input_dir"/*$forward_file_suffix; do
    metawrap binning "$binner" $forward_file $(echo "$forward_file" | sed "s/$forward_file_suffix/$reverse_file_suffix/") \
     --single-end $(echo "$forward_file" | sed "s/$forward_file_suffix/$forward_unpaired_file_suffix/") \
    $(echo "$forward_file" | sed "s/$forward_file_suffix/$reverse_unpaired_file_suffix/" --interleaved \
    -a "$assembly" -o $(echo "$binner" | sed 's/--//g') -t "$threads" -m "$ram" )
  done
done
