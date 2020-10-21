#!/bin/bash
# Pipeline installation script
# Written by: Phagomica group
# Last updated on: 2020-09-06

# Create Conda environments
echo "Creating hummann2 environment..."
conda env create --name humann2 --file ../conda-envs/humann2-env.yml

echo "Creating preprocessing environment..."
conda env create --name preprocessing --file ../conda-envs/preprocessing.yml

# TODO crear ambientes de conda y archivos yml

#conda env create --name binning --file # TODO poner yml

#conda env create --name assembly --file # TODO poner yml

#conda env create --name metaphlan --file # TODO poner yml

#conda env create --name picking16S --file # TODO poner yml
