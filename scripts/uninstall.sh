#!/bin/bash
# Pipeline uninstallation script
# Written by: Phagomica group
# Last updated on: 2020-10-22

SCRIPTS_DIR=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
source "$SCRIPTS_DIR"/config.sh

# Remove Conda environments
remove_envs(){
  echo "Removing preprocessing environment..."
  conda env remove --name metabiome-preprocessing

  echo "Removing binning environment..."
  conda env remove --name metabiome-taxonomic-binning

  echo "Removing assembly environment..."
  conda env remove --name metabiome-genome-assembly

  echo "Removing metaphlan environment..."
  conda env remove --name metabiome-taxonomic-profiling

  echo "Removing picking16S environment..."
  conda env remove --name metabiome-picking16S

  echo "Removing metabat2 environment..."
  conda env remove --name metabiome-metabat2

  echo "Removing maxbin2 environment..."
  conda env remove --name metabiome-maxbin2

  echo "Removing concoct environment..."
  conda env remove --name metabiome-concoct
}

remove_links(){
  echo "Uninstalling metabiome..."
  unlink "$HOME"/.local/bin/metabiome
  echo "Unlinking metabiome from $COMPLETION_DIR"
  unlink "$AT_INSTALL_COMPLETION_DIR"/metabiome
}

while true; do
  read -p "Are you sure to uninstall Metabiome? [y/n] " answer
  case "$answer" in
    [Yy]* ) remove_envs; remove_links
            sed -i.bak '/AT_INSTALL_COMPLETION_DIR/d' "$SCRIPTS_DIR"/config.sh
            rm "$SCRIPTS_DIR"/config.sh.bak; break;;
    [Nn]* ) exit;;
    * )     echo "Please answer yes or no.";;
  esac
done
