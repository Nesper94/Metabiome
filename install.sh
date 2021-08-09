#!/bin/bash
# Pipeline installation script
# Written by: Phagomica group
# Last updated on: 2020-10-22

METABIOME_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$METABIOME_DIR"/scripts/config.sh

# Create link to metabiome.sh
create_link(){
    if [[ ! -d "$HOME"/.local/bin/ ]]; then
        mkdir -p ~/.local/bin/
    fi

    # Check if .local/bin is in PATH
    echo "$PATH" | grep -q "$HOME/.local/bin"
    if [[ "$?" -ne 0 ]]; then
        echo "PATH=$PATH:$HOME/.local/bin" >> "$HOME"/.bash_profile
    fi

    echo "Creating symlink in $HOME/.local/bin/"
    ln -s "$METABIOME_DIR"/scripts/metabiome.sh "$HOME"/.local/bin/metabiome
}

# Create Conda environments
create_envs(){
    echo "Creating preprocessing environment..."
    conda env create --file "$METABIOME_DIR"/conda-envs/preprocessing.yml

    echo "Creating read-based binning environment..."
    conda env create --file "$METABIOME_DIR"/conda-envs/taxonomic-binning.yml

    echo "Creating assembly environment..."
    conda env create --file "$METABIOME_DIR"/conda-envs/genome-assembly.yml

    echo "Creating profiling environment..."
    conda env create --file "$METABIOME_DIR"/conda-envs/taxonomic-profiling.yml

    echo "Creating extract16S environment..."
    conda env create --file "$METABIOME_DIR"/conda-envs/extract16S.yml

    echo "Creating MetaBAT2 environment..."
    conda env create --file "$METABIOME_DIR"/conda-envs/metabat2.yml

    echo "Creating MaxBin2 environment..."
    conda env create --file "$METABIOME_DIR"/conda-envs/maxbin2.yml

    echo "Creating CONCOCT environment..."
    conda env create --file "$METABIOME_DIR"/conda-envs/concoct.yml

    echo "Creating DAS Tool environment..."
    conda env create --file "$METABIOME_DIR"/conda-envs/das_tool.yml

    echo "Creating CheckM environment..."
    conda env create --file "$METABIOME_DIR"/conda-envs/checkm.yml
}

# Copy link to bash completion script
install_completion(){
    if [[ ! -d "$COMPLETION_DIR" ]]; then
        mkdir -p "$COMPLETION_DIR"
    fi
    echo "Installing completion script..."
    ln -s "$METABIOME_DIR"/scripts/_metabiome "$COMPLETION_DIR"/metabiome
    echo "AT_INSTALL_COMPLETION_DIR=$COMPLETION_DIR" >> "$METABIOME_DIR"/scripts/config.sh
}

create_envs
create_link
install_completion
