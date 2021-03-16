.. _tutorial:

Tutorial
========

The purpose of this tutorial is to perform several steps of a metagenomic
analysis using a data set from `Project PRJEB10295 <https://www.ebi.ac.uk/ena/browser/view/PRJEB10295>`_, through our pipeline
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

The  data set for this tutorial is from the project *PRJEB10295*, which is
a metagenomic study of the human palms. It consists of two samples derived
from paired-end sequencing: *ERR981212* and *EEE981213*. Without further ado,
let's download these metagenomic samples like so:

.. code-block:: bash

    # Create directory of the raw reads
    mkdir sample_data
    # Download raw reads for downstream analysis
    wget -P sample_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR981/ERR981212/ERR981212_1.fastq.gz
    wget -P sample_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR981/ERR981212/ERR981212_2.fastq.gz
    wget -P sample_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR981/ERR981213/ERR981213_1.fastq.gz
    wget -P sample_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR981/ERR981213/ERR981213_2.fastq.gz

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
:code:`bowtie2` command. But before running :code:`bowtie2`, we will need to 
subsample the Human reference genome, for tutorial purposes only: 

.. code-block:: bash

    # Activate environment to subsample Human Reference Genome
    conda activate metabiome-preprocessing
    # Download Human Reference Genome
    wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
    # Decompress Human Reference Genome
    gunzip GRCh38_latest_genomic.fna.gz
    # Subsample Human Reference Genome
    fasta-subsample GRCh38_latest_genomic.fna 1 -norand > GRCh38_sub.fna

.. note:: Be aware that we subsampled the Human Reference Genome in order to 
    perform the decontamination step quickly and smoothly. However, for real 
    metagenomic studies you should always use the whole Human Reference Genome.

Now that we have subsampled the Human Reference Genome, let's perform the decontamination with :code:`bowtie2` command like so:

.. code-block:: bash

    metabiome bowtie2 -i filtered_reads/ -o decontaminated_reads/ -hu GRCh38_sub.fna 

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

Using Kaiju
...........

First, let's do it through :code:`kaiju` command. This command will perform the
taxonomic binning, but focusing only in viral communities from your metagenomic
samples.

.. code-block:: bash

    metabiome kaiju -i decontaminated_reads/ -o kaiju_out/ -x taxa_names/ -k krona/ -D kaiju_db/ -d viruses

From this running, you will find two main output directories:
:file:`taxa_names/` and :file:`krona/`, which contain the taxa classification of
the assigned reads and their visualization through krona figures, respectively.

Using Kraken
............

To perform the taxonomic binning with Kraken, we must first download a database
for Kraken to use. In `this link <https://benlangmead.github.io/aws-indexes/k2>`_
you can find a set of different databases to use with Kraken depending on your
needs. In this tutorial, we will use the Viral database just because it is a
lightweight one and you can download it quickly:

.. code-block:: bash

    # Download and extract Viral database
    mkdir kraken2_db
    cd kraken2_db
    wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20201202.tar.gz
    tar -xvzf k2_viral_20201202.tar.gz

Now that we have a database, we can perform the taxonomic classification using
the following command:

.. code-block:: bash

    ../scripts/kraken2.sh -i decontaminated_reads/ -o kraken2_out/ -db kraken2_db/

Visualizing Kraken results
''''''''''''''''''''''''''

We have just performed the taxonomic classification of our reads with Kraken, so
let's visualize these results using Krona:

.. code-block:: bash

    metabiome krona -i kraken2_out/ -o krona_out/

And that's all! Inside the :file:`krona_out/` folder you will now find the Krona
graphs displaying the composition of your samples.

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


Using MetaSPAdes
................

.. code-block:: bash

    # metaSPAdes
    metabiome metaspades -i decontaminated_reads/ -o metaspades-assembled-reads/


Using MEGAHIT
.............

.. code-block:: bash

    # MEGAHIT
    metabiome megahit -i decontaminated_reads/ -o megahit-assembled-reads/

These output genome draft assemblies are frequently used to perform genome quality assessment
and binning.

Quality assembly
----------------

Genome binning
**************

The following step is to generate bins from the previous draft genomes or
contigs. To do so, we will use three different binners::code:`Metabat2`,
:code:`Maxbin2` and :code:`CONCOCT`. Let's begin with :code:`Metabat2`, but
before that let's generate a read coverage table with the next command:

Using Metabat2
--------------

.. code-block:: bash
    
    # Generate read coverage table for Metabat2 running
    metabiome coverage_table  -i contigs_reads/ -o read_coverage/

Now, let's run :code:`Metabat2` through the next command:

.. code-block:: bash

    # Metabat2
    metabiome metabat2 -i contigs/ -co read_coverage/ -o metabat2/ 

Using Maxbin2
-------------

The next binner will be :code:`Maxbin2`. Let's run the command like so: 

.. code-block:: bash

    # Maxbin2
    metabiome maxbin2 -i contigs_reads/ -o maxbin2_out/

Using CONCOCT
-------------

Last but not least, let's run :code:`CONCOCT` command:

.. code-block:: bash

    # CONCOCT
    metabiome concoct -i contigs_reads/ -o concoct_out/
