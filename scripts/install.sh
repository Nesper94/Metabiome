#!/bin/bash
# Pipeline installation script
# Written by: Phagomica group
# Last updated on: 2020-10-22

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
