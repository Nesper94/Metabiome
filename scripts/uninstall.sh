#!/bin/bash
# Pipeline uninstallation script
# Written by: Phagomica group
# Last updated on: 2020-10-22

# Remove Conda environments
remove_envs(){
  echo "Removing hummann2 environment..."
  conda env remove --name metabiome-humann2

  echo "Removing preprocessing environment..."
  conda env remove --name metabiome-preprocessing

  echo "Removing binning environment..."
  conda env remove --name metabiome-taxonomic-binning

  echo "Removing assembly environment..."
  conda env remove --name metabiome-genome-assembly

  echo "Removing metaphlan environment..."
  conda env remove --name metabiome-metaphlan

  echo "Removing picking16S environment..."
  conda env remove --name metabiome-picking16S
}

remove_links(){
  echo "Uninstalling metabiome..."
  unlink "$HOME"/.local/bin/metabiome
}

restore_bash_profile(){
  echo "Removing changes in ~/.bash_profile ..."
  sed -i.bak '/>>Metabiome>>/,/<<Metabiome<</d' ~/.bash_profile
  rm ~/.bash_profile.bak
}

while true; do
  read -p "Are you sure to uninstall Metabiome? [y/n] " answer
  case "$answer" in
    [Yy]* ) remove_envs; remove_links; restore_bash_profile; break;;
    [Nn]* ) exit;;
    * )     echo "Please answer yes or no.";;
  esac
done
