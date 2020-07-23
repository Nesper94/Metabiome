#!/bin/bash

# We will follow this tutorial for data processing:
# https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-Tutorial-(Humann2)
# Use FastQC to evaluate quality of reads

sudo apt install fastqc
fastqc ../raw-data/BB190_L001_R1_001.fastq

# Trimm reads with Trimmomatic
# Trimmomatic should be downloaded from the official web: http://www.usadellab.org/cms/?page=trimmomatic
# for using older versions can be troubling.

wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip

cd ../

#This worked better:

read -p "Forward file: " f-file
read -p "Reverse file: " r-file

java -jar trimmomatic-0.39.jar PE raw-data/BB191_L001_R1_001.fastq raw-data/BB191_L001_R2_001.fastq \
trimmed-reads/BB191_f_paired.fq.gz trimmed-reads/BB191_f_unpaired.fq.gz \
trimmed-reads/BB191_r_paired.fq.gz trimmed-reads/BB191_r_unpaired.fq.gz \
TRAILING:30 MINLEN:130 ILLUMINACLIP:nextera.fa:2:30:10:2

# Bash solution to run Trimmomatic in a great number of samples found at
# https://www.biostars.org/p/294842/

for i in `ls -1 *R1*.fastq | sed 's/\_R1.fastq//'`; do echo trimmomatic PE -phred33 $i\_R1.fastq $i\_R2.fastq $i\_R1_paired.fq.gz $i\_R1_unpaired.fq.gz $i\_R2_paired.fq.gz $i\_R2_unpaired.fq.gz ILLUMINACLIP:contams_forward_rev.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >> cmd_file; done

# This will create a file named cmd_file with the commands to run trimmomatic
