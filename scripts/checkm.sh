#!/bin/bash
# Wrapper script to assess bin quality with CheckM.
# Written by: Phagomica Group
# Last updated on: 2021-07-21

set -e
SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/functions.sh

function usage() {
cat <<HELP_USAGE
Assess bin quality with CheckM workflows.
Usage: metabiome checkm [options] -i <in_dir> -o <out_dir> [-opts CheckM_options]

Required:
-i in_dir           Input directory containing bins per samples. (each sample with their
                    respective bins must have its own directory)
-o out_dir          Directory in which results will be saved. This directory
                    will be created if it doesn't exist.

Options:
-w STRING           CheckM's Workflow. Choose between lineage_wf (lineage analysis workflow)
                    or taxonomy_wf (taxonomic analysis workflow). (default=lineage_wf)
-r STRING           Taxonomic rank (Required only if taxonomic_wf is chosen in -w flag)
-ta STRING          Taxon of interest (Required only if taxonomic_wf is chosen in -w flag)
-t  NUM             Number of threads to use. (default=1)
-opts OPTIONS       CheckM's options.
-h, --help          Show this help.
-hh                 Show CheckM's help message.

HELP_USAGE
}

# Exit if command is called with no arguments
validate_arguments "$#"

# Parse command line arguments
while (("$#")); do
    case "$1" in
        -h|--help ) usage; exit 0 ;;
        -hh )      activate_env metabiome-checkm; checkm -h; exit 0 ;;
        -i )       input_dir=$(readlink -f "$2"); shift 2 ;;
        -o )       out_dir=$(readlink -m "$2"); shift 2 ;;
        -w )       workflow="$2"; shift 2 ;;
        -r )       rank="$2"; shift 2 ;;
        -ta )      taxon="$2"; shift 2 ;;
        -t )       threads="$2"; shift 2 ;;
        -opts )    shift; checkm_opts="$@"; break ;;
        * )        echo "Option '$1' not recognized"; exit 1 ;;
    esac
done

# Verify that input directory exists
validate_input_dir

# Create output directory if it doesn't exists
validate_output_dir

# Activate conda environment
activate_env metabiome-checkm

# Output info
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Input directory: ${input_dir:?'Input directory not set'}"
echo "Number of threads: ${threads:=1}"
echo "CheckM workflow: ${workflow:=lineage_wf}"
echo "CheckM called with options: $checkm_opts"

# Run lineage workflow
if [[ "$workflow" == "lineage_wf" ]]; then
    echo "Perform lineage analysis workflow"
    create_dir "$out_dir" "$workflow"

    for sample in "$input_dir"/*; do
        core_name=$(get_core_name "$sample")
        create_dir "$out_dir"/"$workflow" "$core_name"

        # Run CheckM
        checkm "$workflow" \
            "$sample" \
            "$out_dir"/"$workflow"/"$core_name" \
            -t "$threads" \
            -f "$out_dir"/"$workflow"/"${core_name}"_results_summary.txt \
            $checkm_opts
    done

# Run taxonomy workflow
elif [[ "$workflow" == "taxonomy_wf" ]] && [[ -n "$rank" ]] && [[ -n "$taxon" ]]; then
    echo "Perform taxonomic analysis workflow"
    create_dir "$out_dir" "$workflow"

    for sample in "$input_dir"/*; do
        core_name=$(get_core_name "$sample")
        create_dir "$out_dir"/"$workflow" "$core_name"

        # Run CheckM
        checkm "$workflow" \
            "$rank" \
            "$taxon" \
            "$sample" \
            "$out_dir"/"$workflow"/"$core_name" \
            -t "$threads" \
            -f "$out_dir"/"$workflow"/"${core_name}"_results_summary.txt \
            $checkm_opts
    done

elif [[ ! -n "$rank" ]] || [[ ! -n "$taxon" ]]; then
    echo "Rank's (-r) or taxon's (-x) flags were not set"
fi
