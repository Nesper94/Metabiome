.. _tutorial:

Tutorial
========

.. contents::

Getting help
------------

One of the most useful things you should learn is how to get help from
Metabiome. Fortunately, this is quite easy; if you want to get help about
Metabiome itself and its modules just execute:

.. code-block:: bash

    metabiome -h
    # Or
    metabiome --help

If you want to get help about a particular module, for example the :code:`qc`
module, execute:

.. code-block:: bash

    metabiome qc -h
    # Or
    metabiome qc --help

These commands will show you how to use Metabiome and its modules and which
parameters it needs or accepts.

Downloading sample data
-----------------------

.. code-block:: bash

    mkdir sample_data
    wget -P sample_data

Quality check
-------------

Now that we have the data, we are going to check the quality of the reads by
using the command :code:`qc` from Metabiome:

.. code-block:: bash

    metabiome qc -i sample_data -o quality_check

After running this command the folder :file:`quality-check/` will be created
and inside it you will find a FastQC report with quality info for each input
file. You can also view this info summarized in the file from MultiQC.

Quality filtering
-----------------

The info from the quality check can now be used to trim and remove bad quality
positions and reads by using the :code:`trimmomatic` command. In this case we
will keep only reads whose minimum length is 150 base pairs (bp) and then we
will remove the last 20 bp because these have lower quality:

.. code-block:: bash

    metabiome trimmomatic -i sample_data -o filtered_reads -opts MINLEN:150 TRAILING:20

Decontamination
---------------

The next step is to remove contaminant reads from our data. Two common
contaminants are sequences coming from researchers or people manipulating the
samples and sequences from the Phi-X174 phage used as control in the
sequencing machines, so we will remove reads coming from these sources using
:code:`bowtie2`.

.. code-block:: bash

    metabiome bowtie2 -i filtered_reads -o decontaminated_reads -hu -ph

Taxonomic profiling
-------------------

Functional profiling
--------------------

The first time you use HUMAnN, you must download two databases, ChocoPhlAn and
a translated search database (UniRef), see `HUMAnN documentation
<https://github.com/biobakery/humann#5-download-the-databases>`_ for more info
about this. Here we will download the demo version of ChocoPhlAn database and
the demo version of UniRef90 database by running the following commands:

.. code-block:: bash

    # Activate environment containing HUMAnN
    conda activate metabiome-taxonomic-profiling

    # Create folder in which databases will be saved
    mkdir humann-db

    # Download databases
    humann_databases --download chocophlan DEMO humann-db/
    humann_databases --download uniref DEMO_diamond humann-db/

After downloading databases we are ready to profile our samples with HUMAnN:

.. code-block:: bash

    metabiome humann -i decontaminated-reads -o humann-results

Taxonomic binning
-----------------

Assembly
--------

In this step you can use two different assemblers that receive the output from :code:`bowtie2`:
metaSPAdes and MEGAHIT, in order to obtain longer sequences. For this, run the following commands:


.. code-block:: bash

    # metaSPAdes
    metabiome metaspades -i decontaminated-reads -o metaspades-assembled-reads

.. code-block:: bash

    # MEGAHIT
    metabiome megahit -i decontaminated-reads -o megahit-assembled-reads

This resulted sequences are frequently used to know the taxonomic profiling.

Conting binning
---------------

Quality assembly
----------------
