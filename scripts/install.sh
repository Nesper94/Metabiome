#!/bin/bash
# Pipeline installation script
# Written by: Phagomica group
# Last updated on: 2020-09-06

# Create Conda environments
echo "Creating hummann2 environment..."
conda create --name humann2 --file ../conda-envs/humann2-env.yml

echo "Creating preprocessing environment..."
conda create --name preprocessing --file preprocessing.yml

# TODO crear ambientes de conda y archivos yml

#conda create --name binning --file # TODO poner yml

#conda create --name assembly --file # TODO poner yml

#conda create --name metaphlan --file # TODO poner yml

#conda create --name picking16S --file # TODO poner yml
