# Metabiome
A flexible and modular pipeline for microbiome analysis.

[![License: MIT](https://img.shields.io/badge/License-MIT-orange.svg)](https://github.com/Nesper94/Metabiome/blob/master/LICENSE)
[![Maintained?](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/Nesper94/Metabiome/pulse)
[![Made with: bash](https://img.shields.io/badge/Made%20with-Bash-1f425f.svg)](https://www.gnu.org/software/bash/)
[![Read: the docs](https://img.shields.io/badge/read-the%20docs-blue)](https://metabiome.readthedocs.io/en/latest/)

## About Metabiome

`Metabiome` is a bioinformatic pipeline for microbiome analysis. It consists of
a series of wrapper scripts written in Bash. It is also contained in Conda
environments to ease its use. `Metabiome` includes several features of common
microbiome analysis:

### Reads preprocessing
- FastQC
- MultiQC
- Trimmomatic
- Bowtie2

### Taxonomic and functional profiling
- HUMAnN 3.0
- MetaPhlAn 3.0

### Taxonomic binning
- Kraken2
- Kaiju

### Genome assembly
- MetaSPADES
- MEGAHIT

#### Assembly evaluation
- MetaQUAST

### Genome binning
- CONCOCT
- MetaBAT2
- MaxBin2
- CoverM

#### Bins refinement and assessment
- DAS Tool
- CheckM

## Pipeline
![Pipeline](https://i.imgur.com/ZpCIXYV.png)
)

## Prerequisites for running Metabiome

There are two prerequisites for running `Metabiome`:

1. [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed
in Linux systems.
2. Metagenomic samples to be analyzed.

## Installation

Download or clone this repository and run `install.sh` in the command line:

```bash
cd metabiome/
bash install.sh
```
The installation script will create Conda environments with all the software you
can use with `Metabiome`.


## Basic usage

`Metabiome` can be invoked like so:

```bash
metabiome --help
```
See the  [documentation](https://metabiome.readthedocs.io/en/latest/) for
further details.


## Contribute to Metabiome

If you want to contribute to `Metabiome`, you should fork the project.
Then, make your desired changes and submit a pull request in order to be
evaluated by the current developers.

## Submit an issue

If you have problems or suggestions regarding `Metabiome`, please submit an
issue on [Metabiome issues](https://github.com/Nesper94/Metabiome/issues).
Also, please search the previous issues (if any) before submiting a new one in
order to avoid duplications.

## Ask for help

If you have questions about how to use `Metabiome`, or are interested in
active collaboration of this project, do not hesitate to contact:

Cristian Grisales-Vargas: cristian.grisales@udea.edu.co
Estefany Lorenzana: estefany.lorenzana@udea.edu.co
Juan Camilo Arboleda: juan.arboleda2@udea.edu.co

## Acknowledgements

Authors: Cristian Grisales-Vargas, Estefany Lorenzana, Juan Camilo Arboleda,
Juan Esteban Pérez Jaramillo.
