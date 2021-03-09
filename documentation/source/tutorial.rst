.. _tutorial:

Tutorial
========

The purpose of this tutorial is to perform several steps of a metagenomic
analysis using a data set from (), through our pipeline
`Metabiome <https://github.com/Nesper94/metabiome>`_ .

By the end of the tutorial, you will be able to:
    * Get to know the ``Metabiome`` working environment.
    * Check the quality of metagenomic reads.
    * Filter and decontaminate metagenomic reads.
    * Perform the taxonomic profiling of metagenomic reads.
    * Perform the taxonomic binning of metagenomic reads.
    * Perform the functional profiling of metagenomic reads.
    * Pick 16S rDNA sequences from metagenomic reads.
    * Assembly metagenomic paired-end reads into contigs.
    * Assess the quality of the metagenomic contigs.
    * Generate bins with metagenomic contigs and their respective paired-end reads.

.. contents::

Getting help
************

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

Tutorial Data Set
*****************

.. code-block:: bash

    mkdir sample_data
    wget -P sample_data

Preprocessing
*************

Quality check
-------------

Now that we have the data, we are going to check the quality of the reads by
using the command :code:`qc` from Metabiome:

.. code-block:: bash

    metabiome qc -i sample_data/ -o quality_check/

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

    metabiome trimmomatic -i sample_data/ -o filtered_reads/ -opts MINLEN:150 TRAILING:20

Decontamination
---------------

The next step is to remove contaminant reads from our data. Two common
contaminants are sequences coming from researchers or people manipulating the
samples and sequences from the Phi-X174 phage used as control in the
sequencing machines, so we will remove reads coming from these sources using
:code:`bowtie2` command. To obtain these sequences, you must download the Human
Reference Genome (`GRCh38.p13 <https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39>`_)
and the PhiX (`phi-X174 <https://www.ncbi.nlm.nih.gov/nuccore/9626372>`_) Genome
from the provided links.

Now that we have downloaded the human and phage reference genomes,
let's perform the decontamination with :code:`bowtie2` command like so:

.. code-block:: bash

    metabiome bowtie2 -i filtered_reads/ -o decontaminated_reads/ -hu Human.fasta -ph PhiX_NC_001422.1.fasta

The most important output files from this step are located in
:file:`decontaminated_reads/`. These files are each of the paired-end and
single-end reads in gzip format, and the summary stats from the alignments.
For example, assume your output file prefix is output:

+-------------------------------------+--------------------------------------------------------------+
| File                                | Description                                                  |
+=====================================+==============================================================+
| (output)_paired_bt2_1.fq.gz         | decontaminated forward paired-end reads in gzipped format.   |
+-------------------------------------+--------------------------------------------------------------+
| (output)_paired_bt2_2.fq.gz         | decontaminated reverse paired-end reads in gzipped format.   |
+-------------------------------------+--------------------------------------------------------------+
| (output)_paired_bt2_summary.txt     | summary stats for paired-end alignment.                      |
+-------------------------------------+--------------------------------------------------------------+
| (output)_unpaired_bt2_f.fq.gz       | Decontaminated forward single-end reads in gzipped format.   |
+-------------------------------------+--------------------------------------------------------------+
| (output)_unpaired_bt2_f_summary.txt | summary stats for forward single-end alignment.              |
+-------------------------------------+--------------------------------------------------------------+
| (output)_unpaired_bt2_r.fq.gz       | Decontaminated reverse single-end reads in gzipped format.   |
+-------------------------------------+--------------------------------------------------------------+
| (output)_unpaired_bt2_r_summary.txt | summary stats for reverse single-end alignment.              |
+-------------------------------------+--------------------------------------------------------------+

Read-based analysis
*******************

Taxonomic profiling
-------------------

Now, consider that you want to predict the taxonomic identity and relative
abundance of your metagenomic samples. To do so, run the :code:`metaphlan3`
command like so:

.. code-block:: bash

    metabiome metaphlan3 -i decontaminated_reads/ -o mphlan_out/

In the ouput directory :file:`mphlan_out/`, you will find the taxa identity and
relative abundances from the metagenomic samples.


Taxonomic binning
-----------------

In addition to taxonomic profiling, you can also predict the taxonomic identity
of your metagenomic samples by taxonomic binning. You can perform the taxonomic
binning through :code:`kaiju` or :code:`kraken2` commands.

First, let's do it through :code:`kaiju` command. This command will perform the
taxonomic binning, but focusing only in viral communities from your metagenomic
samples.

.. code-block:: bash

    metabiome kaiju -i decontaminated_reads/ -o kaiju_out/ -x taxa_names/ -k krona/ -D kaiju_db/ -d viruses

From this running, you will find two main output directories:
:file:`taxa_names/` and :file:`krona/`, which contain the taxa classification of
the assigned reads and their visualization through krona figures, respectively.


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
    mkdir humann_db

    # Download databases
    humann_databases --download chocophlan DEMO humann_db/
    humann_databases --download uniref DEMO_diamond humann_db/

After downloading databases we are ready to profile our samples with HUMAnN:

.. code-block:: bash

    metabiome humann -i decontaminated_reads/ -o humann-results/


16S rDNA picking
----------------
Now, lets suppose you want to perform additional analyses based on the 16S rDNA.
The :code:`BBDuk` command can pick the 16S rDNA from your metagenomic samples.
But first, you will need to download the 16S rDNA sequences from the database of
your choice. We recommend to download the 16S rDNA sequences from the up-to-date
`SILVA_16S database <https://www.arb-silva.de/>`_ and store it in a directory
(:file:`SILVA_16S/`)

.. code-block:: bash

    metabiome bbduk -i decontaminated_reads/ -o bbduk_out/ -D SILVA_16S/

The output of :code:`BBDuk` is located in :file:`bbduk_out/`. This output is
very similar to the `Decontamination section <Decontamination_>`_ output.
However, in this context, these files represent the metagenomic reads that did
aligned to the 16S rDNA sequences.

*De-novo* Assembly
******************

Genome assembly
---------------

In this step you can use two different assemblers that receive the output from
:code:`bowtie2`: metaSPAdes and MEGAHIT, in order to obtain longer sequences.
You can use just the assembler you like the most, or use both as we will do in
this tutorial. To perform the assembly, just run the following commands:

.. code-block:: bash

    # metaSPAdes
    metabiome metaspades -i decontaminated_reads/ -o metaspades-assembled-reads/

    # MEGAHIT
    metabiome megahit -i decontaminated_reads/ -o megahit-assembled-reads/

This resulted sequences are frequently used to know the taxonomic profiling.

Quality assembly
----------------

Genome binning
**************

The following step is to generate bins from the previous draft genomes or
contigs. To do so, we will use three different binners::code:`Metabat2`,
:code:`Maxbin2` and :code:`CONCOCT`. Let's begin with :code:`Metabat2`, but
before that let's generate a read coverage table with the next command:

.. code-block:: bash
    
    ##Generate read coverage table for Metabat2 running:
    metabiome coverage_table  -i contigs_reads/ -o read_coverage/

Now, let's run :code:`Metabat2` through the next command:

.. code-block:: bash

    ##Metabat2
    metabiome metabat2 -i contigs/ -co read_coverage/ -o metabat2/ 

The next binner will be :code:`Maxbin2`. Let's run the command like so: 

.. code-block:: bash

    ##Maxbin2
    metabiome maxbin2 -i contigs_reads/ -o maxbin2_out/

Last but not least, let's run :code:`CONCOCT` command:

.. code-block:: bash

    ##CONCOCT
    metabiome concoct -i contigs_reads/ -o concoct_out/
