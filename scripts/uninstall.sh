#!/bin/bash
# Pipeline uninstallation script
# Written by: Phagomica group
# Last updated on: 2020-10-22

# Remove Conda environments
echo "Removing hummann2 environment..."
conda env remove --name hummann2

echo "Removing preprocessing environment..."
conda env remove --name preprocessing

echo "Removing binning environment..."
conda env remove --name binning

echo "Removing assembly environment..."
conda env remove --name assembly

echo "Removing metaphlan environment..."
conda env remove --name metaphlan

echo "Removing picking16S environment..."
conda env remove --name picking16S
