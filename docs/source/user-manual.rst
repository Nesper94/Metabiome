.. _usermanual:

User manual
===========

.. contents::

Basic requirements
-------------------

Metabiome was built under different bioinformatics considerations to simplify its use. Prior to installation you require:


1. Linux terminal. If you are a windows user, you can use a virtual machine. For more information, please check `VirtualBox <https://www.virtualbox.org/>`_.
2. CONDA as a package manager. We recommend installing MINICONDA as it contains all necessary packages, and the latest version can be downloaded from `this page <https://docs.conda.io/en/latest/miniconda.html#linux-installers>`_. To install it, run the following commands:

.. code-block:: bash

    bash Miniconda3-latest-Linux-x86_64.sh

3. Access to a computer server, given the large size of metagenomic data sets and the complexity of subsequent analyzes.

Installation
------------

In order to install, clone Metabiome into your machine and execute
:file:`install.sh` by running the following commands:

.. code-block:: bash

    git clone -b master https://github.com/Nesper94/metabiome.git
    cd Metabiome/
    bash install.sh

You can also download Metabiome directly from the `Github page
<https://github.com/Nesper94/metabiome>`_.

The installation script will do the following:

- Create a link to the Metabiome's main script in :file:`~/.local/bin/` and if
  necessary add this folder to :envvar:`PATH` so that you can use Metabiome
  from whatever location in your terminal.
- Create Conda environments with all the software you can use with Metabiome.
- Install the :file:`_metabiome` bash completion script to make Metabiome's
  commands easier to use.

Uninstalling
------------

If for some reason you need to uninstall Metabiome, all you have to do is go to
the Metabiome's folder and execute the :file:`uninstall.sh` script:

.. code-block:: bash

    cd scripts/
    bash uninstall.sh

.. _opts-flag:

Additional command line options
-------------------------------

Some Metabiome commands accept the flag ``-opts``, if this flag is used
then it must be the last flag used in the command. After this flag you can
write options from the native program executed by Metabiome. For example, if
you want to assemble a large and complex metagenome like soil and need to use
the flag ``--presets`` from MEGAHIT, then you would do the following:

.. code-block:: bash

    metabiome megahit -i in_dir -o out_dir -opts --presets meta-large


Naming convention
-----------------

Metabiome's working needs that the input files agree with one of the two most
widespread naming conventions for paired end sequences:
Illumina :file:`_R1_`/:file:`_R2_` naming convention and :file:`_1`/:file:`_2`
naming convention.

File extensions accepted are:

- :file:`.fq`
- :file:`.fastq`
- :file:`.fq.gz`
- :file:`.fastq.gz`

.. _modules:

Metabiome modules
-----------------

Metabiome contains 9 modules that comprise the necessary tools for the analysis
of the main points within metagenomics. They are separated by conda environments,
created from a ``.yaml`` file, which describes the software that each one implements
and the required version. These files are stored in the :file:`conda_envs/` directory. A module
can have one or more than one tool or software, and each one has a separate script
to execute it. They are stored in the :file:`scripts/` directory.

The following is the list of software contained in each one of the Metabiome
modules:

metabiome-preprocessing
***********************

- Bowtie 2 (v2.3)
- FastQC (v0.11)
- MultiQC (v1.6)
- Trimmomatic (v0.39)

metabiome-taxonomic-profiling
*****************************

- MetaPhlAn3 (v3.0)
- HUMAnN3 (v3.0)

metabiome-taxonomic-binning
***************************

- Kaiju (v1.7)
- Kraken2 (v2.1)
- Krona (v2.7)

metabiome-extract16S
********************

- BBMap (v38.87)

metabiome-genome-assembly
*************************

- MEGAHIT (v1.2)
- QUAST (v5.0)
- SPAdes (v3.12)

metabiome-concoct
*****************

- CONCOCT (v1.1)
- CoverM (v0.6)

metabiome-maxbin2
*****************

- MaxBin2 (v2.2)

metabiome-metabat2
******************

- MetaBAT2 (v2.15)

metabiome-das_tool
******************

- DAS Tool (v1.1.2)

Debug mode
----------

Sometimes errors can appear when you try to run a command. If you are a
developer and/or have knowledge in bash, you can use the Metabiome's debug
mode to obtain more information about the error. The debug mode is turned on
by setting the :envvar:`DEBUG_METABIOME` environment variable:

.. code-block:: bash

    export DEBUG_METABIOME="yes"

Then you can run the command again and see a detailed output that can help
in the debugging process.

Software native help message
-----------------------------

Metabiome commmand line option includes an additional flag named ``-hh``
that allows the users to see the programs native help message. For example,
if you want more information about the assembler metaSPAdes the ``-hh``
would show the help message from this software:

.. code-block:: bash

    metabiome metaspades -hh
    

.. note::

    It seems that Krona doesn't have a help command, and FastQC 
    and MultiQC usage in our pipeline is very simple, so they don't have this option.
