#!/bin/bash
# Pipeline installation script
# Written by: Phagomica group
# Last updated on: 2020-10-22

METABIOME_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$METABIOME_DIR"/scripts/config.sh

create_link(){
    echo "Creating symlink in $HOME/.local/bin/"
    ln -s "$METABIOME_DIR"/scripts/metabiome.sh "$HOME"/.local/bin/metabiome
}

# Create link to metabiome.sh
echo "$PATH" | grep -q "$HOME/.local/bin"
if [[ "$?" -eq 0 ]]; then
    if [[ ! -d "$HOME"/.local/bin/ ]]; then
        mkdir -p "$HOME"/.local/bin/
    fi
    create_link
else
    while true; do
        read -p "$HOME/.local/bin is not in your PATH, do you want to add it? [y/n] " answer
        case "$answer" in
            [Yy]* ) if [[ ! -d "$HOME"/.local/bin/ ]]; then
                        mkdir -p "$HOME"/.local/bin/
                    fi
                    echo "PATH=$PATH:$HOME/.local/bin" >> "$HOME"/.bashrc
                    create_link; break;;
            [Nn]* ) echo "Aborting installation..."; exit;;
            * )     echo "Please answer yes or no.";;
        esac
    done
fi

# Create Conda environments
echo "Creating hummann2 environment..."
conda env create --file "$METABIOME_DIR"/conda-envs/humann2-env.yml

echo "Creating preprocessing environment..."
conda env create --file "$METABIOME_DIR"/conda-envs/preprocessing.yml

echo "Creating read-based binning environment..."
conda env create --file "$METABIOME_DIR"/conda-envs/read-binning.yml

echo "Creating assembly environment..."
conda env create --file "$METABIOME_DIR"/conda-envs/assembly.yml

echo "Creating MetaPhlAn3 environment..."
conda env create --file "$METABIOME_DIR"/conda-envs/metaphlan.yml

#conda env create --file # TODO: poner yml picking16S

# Source completion script
echo "Installing completion script..."

cat << EOF >> ~/.bash_profile

# >>Metabiome>>
# Contents in this block are managed by Metabiome
if [[ -f "$METABIOME_DIR"/scripts/_metabiome ]]; then
    . "$METABIOME_DIR"/scripts/_metabiome
fi
# <<Metabiome<<
EOF
