.. _tutorial:

Tutorial
========

.. contents::

Downloading sample data
-----------------------

.. code-block:: bash

    mkdir sample-data
    wget -P sample-data

Quality Check
-------------

.. code-block:: bash

    metabiome qc -i sample-data -o quality-check

Decontamination
---------------

Taxonomic profiling
-------------------

Functional profiling
--------------------

The first time you use HUMAnN, you must download two databases, ChocoPhlAn and
a translated search database (UniRef), see `HUMAnN documentation
<https://github.com/biobakery/humann#5-download-the-databases>`_ for more info
about this. Here we will download the full ChocoPhlAn database and the full
UniRef90 database by running the following commands:

.. code-block:: bash

    # Activate environment containing HUMAnN
    conda activate metabiome-taxonomic-profiling

    # Create folder in which databases will be saved
    mkdir humann-db

    # Download databases
    humann_databases --download chocophlan full humann-db/
    humann_databases --download uniref uniref90_diamond humann-db/

After downloading databases we are ready to profile our samples with HUMAnN:

.. code-block:: bash

    metabiome humann -i decontaminated-reads -o humann-results

Taxonomic binning
-----------------

Assembly
--------