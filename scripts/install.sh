#!/bin/bash
# Pipeline installation script
# Written by: Phagomica group
# Last updated on: 2020-10-22

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/config.sh

# Create Conda environments
echo "Creating hummann2 environment..."
conda env create --file ../conda-envs/humann2-env.yml

echo "Creating preprocessing environment..."
conda env create --file ../conda-envs/preprocessing.yml

# TODO crear ambientes de conda y archivos yml

#conda env create --file # TODO poner yml binning

echo "Creating assembly environment..."
conda env create --file ../conda-envs/assembly.yml

#conda env create --file # TODO poner yml metaphlan

#conda env create --file # TODO poner yml picking16S

# Copy link to bash completion script
if [[ ! -d "$COMPLETION_DIR" ]]; then
  mkdir -p "$COMPLETION_DIR"
fi
ln -s "$SCRIPTS_DIR"/_metabiome "$COMPLETION_DIR"/metabiome
