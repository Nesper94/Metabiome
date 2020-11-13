#!/bin/bash
# Trimmomatic wrapper script to make quality control on FASTQ files
# Last updated on: 2020-08-29

set -e

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
    echo "Usage: $0 -i <input directory> -o <output directory> [-t <threads>] ['TRIMMOMATIC_OPTIONS']"
    echo ""
    echo "Make sure TRIMMOMATIC_OPTIONS are enclosed within quotation marks."
    echo "Output directory will be created if it doesn't exists."
}

# Exit if command is called with no arguments
validate_arguments "$#"

while [[ -n "$1" ]]; do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift ;;
        -o )        out_dir=$(readlink -f "$2"); shift ;;
        -t )        threads="$2"; shift ;;
        * )         trimopt="$@" ;;
    esac
    shift
done

# Verify that input directory is set and exists
validate_input_dir

# Create output directory if it doesn't exists.
validate_output_dir

# Activate conda environment
activate_env preprocessing

# Output info
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: ${input_dir}"
echo "Output directory: ${out_dir}"
echo "Number of threads: ${threads:=4}"
echo "Trimmomatic version: $(trimmomatic -version)"

# Las opciones utilizadas con los reads de prueba el 2020-07-13
# fueron: MINLEN:140 TRAILING:25 HEADCROP:20
# Es importante recordar que en Trimmomatic el orden de las opciones indica
# su orden de ejecución.

for file in "$input_dir"/*; do

    # Make sure to process only fastq, fq.gz or fastq.gz files
    if [[ $file == *R1*.fastq ]] || [[ $file == *R1*.fq.gz ]] || [[ $file == *R1*.fastq.gz ]]; then
        trimmomatic PE -threads ${threads:=4} "$file" $(echo "$file" | sed 's/R1/R2/') \
        "$out_dir"/$(basename "$file" | sed 's/R1.*/1_paired_trim.fq.gz/')   \
        "$out_dir"/$(basename "$file" | sed 's/R1.*/1_unpaired_trim.fq.gz/') \
        "$out_dir"/$(basename "$file" | sed 's/R1.*/2_paired_trim.fq.gz/')   \
        "$out_dir"/$(basename "$file" | sed 's/R1.*/2_unpaired_trim.fq.gz/') \
        $trimopt 2>&1 | tee -a "$out_dir"/trimmomatic.log

        echo "" # Add new line in order to make output easier to read

    elif [[ $file == *.fastq ]] || [[ $file == *.fq.gz ]]; then
        :
    else
        echo "$file will not be processed as is not a .fastq or .fq.gz file."
    fi
done

echo "Done. You should now execute Bowtie2 in order to clean contaminant reads."
