# Exit if command is called with no arguments
validate_arguments(){
    if [[ "$1" == 0 ]]; then
        echo "Error: No arguments given." >&2
        usage
        exit 1
    fi
}

# Verify that input directory exists
validate_input_dir(){
    if [[ -z "$input_dir" ]]; then
        echo "$0: Error: Input directory is not set." >&2
        exit 1
    fi

    if [[ ! -d "$input_dir" ]]; then
        echo "$0: Error: $input_dir is not a valid directory." >&2
        exit 1
    fi
}

# Create output directory if it doesn't exists.
validate_output_dir(){
    if [[ -z "$out_dir" ]]; then
        echo "$0: Error: Output directory is not set." >&2
        exit 1
    fi

    if [[ ! -d "$out_dir" ]]; then
        mkdir -p "$out_dir"
    fi
}

# Activate conda environment
activate_env(){
    source activate "$1"
    if [[ "$?" -ne 0 ]]; then
        echo "$0: Error: Make sure Conda base environment is active by running 'conda activate base'" >&2
        exit 1
    fi
}

# Change filename from forward to reverse
forward_to_reverse(){
    if (( $# == 0 )); then
        sed 's/_R1_/_R2_/ ; s/_1\./_2\./' < /dev/stdin
    else
        echo "$1" | sed 's/_R1_/_R2_/ ; s/_1\./_2\./'
    fi
}

remove_forward_suffix(){
    if (( $# == 0 )); then
        sed 's/_R1_/_/ ; s/_1\./\./' < /dev/stdin
    else
        echo "$1" | sed 's/_R1_/_/ ; s/_1\./\./'
    fi
}

get_core_name(){
    if (( $# == 0 )); then
        local filename=$(cat /dev/stdin | xargs basename -- | remove_forward_suffix)
        echo "${filename%%.*}"
    else
        local filename=$(basename "$1" | remove_forward_suffix)
        echo "${filename%%.*}"
    fi
}

create_dir(){
	if (( $# == 2 ));then
        if [[ ! -d "$1"/"$2" ]];then
            echo -e "output directory path: $1"
            echo -e "Name of the folder to be created: $2"
			mkdir "$1"/"$2"
        fi
    fi
}
