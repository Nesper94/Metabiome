#!/bin/bash
# Install SRA Toolkit and data

set -e

sudo apt install sra-toolkit

# Download the accession list from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SAMN09762211&o=acc_s%3Aa
# and save it to raw-data/SRR_Acc_List.txt

# Download the data from the article

cd ../raw-data/

prefetch --option-file SRR_Acc_List.txt

# Data will be downloaded to $HOME/ncbi/public/sra/, then we must convert the
# .sra files to .fastq

fastq-dump ~/ncbi/public/sra/*

# Remove the .sra files

rm $HOME/ncbi/public/sra/*.sra
mv $HOME/ncbi/public/sra/* . # move files to current folder (raw-data/)
