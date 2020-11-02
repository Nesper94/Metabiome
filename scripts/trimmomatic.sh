#!/bin/bash
# Este script instala y ejecuta Trimmomatic usando como input una carpeta con
# archivos en formato FastQ con los reads en bruto.
# Last updated on: 2020-08-29

set -e

function usage() {
    echo "Usage: $0 -i <input directory> -o <output directory> [-t <threads>] 'TRIMMOMATIC_OPTIONS'"
    echo "Make sure TRIMMOMATIC_OPTIONS are enclosed with quotation marks."
    echo "Output directory will be created if it doesn't exists."
}

if [[ "$#" == 0 ]]; then
    echo "Error: No arguments given." >&2
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
        * )         trimopt="$@"
            ;;
    esac
    shift
done

if [ ! -d "$input_dir" ]; then
   echo "Error: You need to specify the directory containing reads in FastQ" >&2
   exit 1
fi

echo "Your read directory is $input_dir"

# Activate conda environment
source activate preprocessing

if [ ! -d "$out_dir" ]; then
    echo "Creating directory '$out_dir'"
	mkdir -p "$out_dir"
fi

# Las opciones utilizadas con los reads de prueba el 2020-07-13
# fueron: MINLEN:140 TRAILING:25 HEADCROP:20
# Es importante recordar que en Trimmomatic el orden de las opciones indica
# su orden de ejecución.

for file in "$input_dir"/*; do
    if [[ $file == *R1*.fastq ]] || [[ $file == *R1*.fq.gz ]]; then  # Make sure to process only fastq or fq.gz files
        trimmomatic PE -threads ${threads:=4} "$file" $(echo "$file" | sed 's/R1/R2/') \
        "$out_dir"/$(basename "$file" | sed 's/R1.*/1_paired_trim.fq.gz/')   \
        "$out_dir"/$(basename "$file" | sed 's/R1.*/1_unpaired_trim.fq.gz/') \
        "$out_dir"/$(basename "$file" | sed 's/R1.*/2_paired_trim.fq.gz/')   \
        "$out_dir"/$(basename "$file" | sed 's/R1.*/2_unpaired_trim.fq.gz/') \
        $trimopt 2>&1 | tee -a "$out_dir"/trimmomatic.log
    elif [[ $file == *.fastq ]] || [[ $file == *.fq.gz ]]; then
        :
    else
        echo "$file will not be processed as is not a .fastq or .fq.gz file."
    fi
done

echo "Done. You should now execute Bowtie2 in order to clean contaminant reads."
