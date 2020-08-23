#!/bin/bash
# Este script instala y ejecuta Trimmomatic usando como input una carpeta con
# archivos en formato FastQ con los reads en bruto.
# Last updated on: 2020-08-21

set -e

function usage() {
    echo "Usage: $0 -i <input directory> -o <output directory> [-t <threads>] [--trim-dir <Trimmomatic directory>]"
    echo "Output directory will be created if it doesn't exists."
}

if [[ "$#" == 0 ]]; then
    echo "No arguments given."
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
        -t )        threads="$2"
            shift
            ;;
        --trim-dir )    trim_loc=$(readlink -f "$2")
            shift
            ;;
        * )         echo "Option '$1' not recognized"; exit 1
            ;;
    esac
    shift
done

script_dir=$(readlink -f $(dirname "$0")) #Lines to get path to main directory
main_dir=$(cd $script_dir ; cd .. ; pwd)

if [ ! -d "$input_dir" ]; then
   echo "Error: You need to specify the directory containing reads in FastQ"
   exit 1
fi

echo "Your read directory is $input_dir"

cd $main_dir

if [ ! -d "$out_dir" ]; then
    echo "Creating directory '$out_dir'"
	mkdir -p "$out_dir"
fi

# Install Trimmomatic

if [[ -z "$trim_loc" ]]; then
    echo 'Installing Trimmomatic...'
    wget -nd -P software/ http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip || exit 1
    cd software/
    unzip Trimmomatic-0.39.zip
    cd $main_dir
fi

trimmomatic=$(ls "$trim_loc"/trimmomatic*.jar)

read -p "Write arguments for Trimmomatic: " trimopt
# Las opciones utilizadas con los reads de prueba el 2020-07-13
# fueron: MINLEN:140 TRAILING:25 HEADCROP:20
# Es importante recordar que en Trimmomatic el orden de las opciones indica
# su orden de ejecución.

cd $input_dir

for i in $(ls *R1*.fastq -1);
do java -jar $trimmomatic PE -threads ${threads:=4} $i $(echo $i | sed 's/R1/R2/') \
"$out_dir"/$(echo $i | sed 's/R1.*/f-paired.fq.gz/')   \
"$out_dir"/$(echo $i | sed 's/R1.*/f-unpaired.fq.gz/') \
"$out_dir"/$(echo $i | sed 's/R1.*/r-paired.fq.gz/')   \
"$out_dir"/$(echo $i | sed 's/R1.*/r-unpaired.fq.gz/') \
$trimopt 2>&1 | tee -a "$out_dir"/trimmomatic.log
done

echo "Done. You should now execute Bowtie2 in order to clean contaminant reads."
