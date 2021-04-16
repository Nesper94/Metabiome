# Metabiome
Pipeline for metagenomics of rhizosphere and associated phages.

## Installation

Download or clone this repository and run `install.sh` in the command line:
```bash
cd metabiome/
bash install.sh
```
The installation script will create Conda environments with all the software you
can use with Metabiome.

## Pre-processing of reads
- FastQC
- MultiQC
- Trimmomatic
- Bowtie2

## Functional profiling
- HUMAnN 2.0 [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/humann2/README.html)
- Kraken2
- Kaiju
- MetaPhlAn3

## Genome assembly
- MetaSPADES
- MEGAHIT

### Evaluation of assembly
- MetaQUAST

## Pipeline
![Pipeline](https://i.imgur.com/nnm7qi3.png)
