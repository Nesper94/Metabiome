.. _usermanual:

User manual
===========

.. contents::

Basic requirements
-------------------

Metabiome was built under different bioinformatics considerations to simplify its use. Prior to installation you require:


1.  Linux terminal. If you are a windows user, you can use a virtual machine. For more information, please check `VirtualBox <https://www.virtualbox.org/>`_.
2. CONDA as a package manager. We recommend installing MINICONDA as it contains all necessary packages, and can be downloaded from `this page <https://docs.conda.io/en/latest/miniconda.html#linux-installers>`_. To install it, run the following commands:

.. code-block:: bash

    sudo apt-get update
    bash Miniconda3-latest-Linux-x86_64.sh

3. Access to a computer server, given the large size of metagenomic data sets and the complexity of subsequent analyzes.

Installation
------------

In order to install, clone Metabiome into your machine and execute
:file:`install.sh` by running the following commands:

.. code-block:: bash

    git clone -b master https://github.com/Nesper94/metabiome.git
    cd metabiome/
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

Use
---

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
