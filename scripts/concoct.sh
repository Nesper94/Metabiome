#!/bin/bash
# CONCOCT wrapper script for the binning of draft assemblies.
# Written by: Phagomica Group
# Last updated on: 2021-02-14

set -e
SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
  cat <<HELP_USAGE
Generate bins from metagenomic samples with CONCOCT.
Usage: metabiome concoct [options] -i <in_dir> -o <out_dir> -opts concoct_options
Required:
-i  in_dir                      Directory containing paired-end reads and contigs in fastq
                                and fasta format, respectively.(The filenames of the contigs and
                                their respective paired-end reads must match)
-o  out_dir                     Directory in which results will be saved. This directory
                                will be created if it doesn't exist.
Options:
-t  NUM                         Number of threads. (default=1)
-ch NUM                         Contigs's chunk size. (default=1000)
-opts concoct_options           CONCOCT's options.
-h, --help                      Show this help.
HELP_USAGE
}
##--------------------------Exiting if input files are missing---------------##:
validate_arguments "$#"
##----------------------Saving input orders into variables------------------##:
while (("$#")); do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -i )        input_dir=$(readlink -f "$2"); shift 2 ;;
        -o )        out_dir=$(readlink -m "$2"); shift 2 ;;
        -ch )       chunk_size="$2"; shift 2 ;;
        -t )        threads="$2"; shift 2 ;;
        -opts )     shift; concoct_opts="$@"; break ;;
        * )         echo "Option '$1' not recognized"; exit 1 ;;
    esac
done
##----------------Verify that input directory exists------------------------##:
validate_input_dir
# Create output directory if it doesn't exists
validate_output_dir
##----------------Activate conda environment--------------------------------##:
activate_env metabiome-concoct
##---------------Output info-------------------------------------------------##:
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Output directory: ${out_dir}"
echo "Number of threads: ${threads:=1}"
echo "CONCOCT version: $(concoct -v)"
echo "Contig chunk size: ${chunk_size:=1000}"
echo "CONCOCT called with options: $concoct_opts"
##-----------------------------CONCOCT binning-------------------------------##:
create_dir "$out_dir" fasta_bins
for file in "$input_dir"/*;do
    #Check correct forward file format:
    if [[ "$file" == *@(*_R1_*|*_1).@(fq|fastq|fq.gz|fastq.gz) ]]; then
            forward_file="$file"
            core_name=$(get_core_name "$forward_file")
            #check the format of the contig or draft genome:
            draft=$(get_genome_format "$input_dir"/"$core_name")
            #create output directory:
            create_dir "$out_dir" "$core_name" && cd "$out_dir"/"$core_name"
            #cut contigs into smaller parts
            cut_up_fasta.py "$draft" -m -o 0 -c "$chunk_size" > "$core_name".k"$chunk_size".fa
            #building Kallisto index:
            [ ! -f "$core_name".idx ] && { echo "Generate Kallisto index:";
            kallisto index "$core_name".k"$chunk_size".fa -i "$core_name".idx;}
            #Map the reads to their contigs from each replicate:
            kallisto quant "$forward_file" $(forward_to_reverse "$forward_file") \
                -i "$core_name".idx -o "$out_dir"/"$core_name" -t "$threads"
            #Generate the coverage table for concoct:
            input_table.py abundance.tsv > "$core_name".kcov.tsv
            #run concoct:
            concoct --composition_file "$core_name".k"$chunk_size".fa \
                --coverage_file "$core_name".kcov.tsv -t "$threads" -b "$core_name" $concoct_opts
            #merge subcontig clustering into original contigs:
            merge_cutup_clustering.py "$core_name"_clustering_gt*.csv > "$core_name"_clust_merged.csv
            create_dir "$out_dir"/fasta_bins "$core_name"
            #extract bins as individual FASTA:
            extract_fasta_bins.py "$draft" "$core_name"_clust_merged.csv \
                --output_path "$out_dir"/fasta_bins/"$core_name"
    fi
done
