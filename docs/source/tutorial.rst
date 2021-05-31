.. _tutorial:

Tutorial
========

The purpose of this tutorial is to perform several steps of a metagenomic
analysis through our pipeline
`Metabiome <https://github.com/Nesper94/metabiome>`_ .

By the end of the tutorial, you will be able to:
    * Get to know the ``Metabiome`` working environment.
    * Check the quality of metagenomic reads.
    * Filter and decontaminate metagenomic reads.
    * Perform the taxonomic profiling of metagenomic reads.
    * Perform the taxonomic binning of metagenomic reads.
    * Perform the functional profiling of metagenomic reads.
    * Extract 16S rDNA sequences from metagenomic reads.
    * Assembly metagenomic paired-end reads into contigs.
    * Assess the quality of the metagenomic contigs.
    * Generate bins with metagenomic contigs and their respective paired-end reads.

.. contents:: Tutorial contents
    :depth: 1
    :local:

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

The  data set for this tutorial is from the
`project PRJEB10295 <https://www.ebi.ac.uk/ena/browser/view/PRJEB10295>`_,
which is a metagenomic study of the human palms. It consists of two samples derived
from paired-end sequencing: *ERR981212* and *EEE981213*. However, for tutorial
purposes only, we have subsampled these files which you can download from here:
`sample data <https://drive.google.com/drive/folders/1TxZPUrRVkoRa8rJNHiOx1sm7GdYN__5y?usp=sharing>`_.
After you download these samples, store them in a directory called
:file:`sample_data` for downstream analysis.

Preprocessing
*************

Quality check
-------------

Now that we have the data, we are going to check the quality of the reads by
using the command :code:`qc` from Metabiome:

.. code-block:: bash

    metabiome qc -i sample_data/ -o quality_check/

After running this command the folder :file:`quality-check/` will be created
and inside it you will find a `FastQC <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_
report with quality info for each input file. You can also view this info
summarized in the file from `MultiQC <https://multiqc.info/>`_.

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
`Bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_. Thus, before
running :code:`bowtie2` command let's download through the next links the
`subsampled Human Genome <https://drive.google.com/file/d/1f49lWDaX63FefH150PZ_p9FUa5UwE5zk/view?usp=sharing>`_
and the `Phi-X174 genome <https://drive.google.com/file/d/1uRdEzysZCySSkBqp-uEn-Cx5MbsQ5F8n/view?usp=sharing>`_,
which we will use to decontaminate the filtered reads like so:


.. warning:: Be aware that we subsampled the Human Reference Genome in order to
    perform the decontamination step quickly and smoothly. However, for real
    metagenomic studies you should always use the whole Human Reference Genome.

.. code-block:: bash

    metabiome bowtie2 -i filtered_reads/ -o decontaminated_reads/ -hu GRCh38_sub.fna \
        -ph PhiX_NC_001422.1.fasta

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
| (output)_unpaired_bt2_f.fq.gz       | decontaminated forward single-end reads in gzipped format.   |
+-------------------------------------+--------------------------------------------------------------+
| (output)_unpaired_bt2_f_summary.txt | summary stats for forward single-end alignment.              |
+-------------------------------------+--------------------------------------------------------------+
| (output)_unpaired_bt2_r.fq.gz       | decontaminated reverse single-end reads in gzipped format.   |
+-------------------------------------+--------------------------------------------------------------+
| (output)_unpaired_bt2_r_summary.txt | summary stats for reverse single-end alignment.              |
+-------------------------------------+--------------------------------------------------------------+


.. warning:: It is important to point out that in this particular case,
    we did not have any reads in the files: :file:`ERR981212_sub_unpaired_bt2_r.fq.gz`
    and :file:`ERR981213_sub_unpaired_bt2_r.fq.gz`. Therefore, we must
    remove these files in order to avoid problems for downstream analysis.
    To do so, take a look at the next command:

    .. code-block:: bash

        rm decontaminated_reads/*sub_unpaired_bt2_r*

Read-based analysis
*******************

Taxonomic profiling
-------------------

Now, consider that you want to predict the taxonomic identity and relative
abundance of your metagenomic samples, through marker-based methods. To do so,
we will use `MetaPhlAn3 <https://huttenhower.sph.harvard.edu/metaphlan/>`_.
However, due to tutorial purposes only, you will have to download our custom
database located here: `metaphlan3_custom_db <https://drive.google.com/drive/folders/1xNzSYTjSYlfycDsSC6_QM47y9Yid9Oe5?usp=sharing>`_.
Be aware that this database is compressed and after downloading it, you must
extract the :file:`metaphlan_custom_db.tar.gz` like so:

.. code-block:: bash

    tar -xvf metaphlan3_custom_db.tar.gz

Now, we can perfom the taxonomic profiling of the metagenomics samples with the
:code:`metaphlan3` command like so:

.. code-block:: bash

    metabiome metaphlan3 -i decontaminated_reads/ -o mphlan_out/ -d metaphlan3_custom_db/ \
        -opts --add_viruses --ignore_eukaryotes --ignore_bacteria --ignore_archaea

In the output directory :file:`mphlan_out/`, you will find the taxa identity and
relative abundances of the metagenomic samples. Additionally, you will find the
following file :file:`merged_mphlan.txt`, which contains the taxonomic profiling
of all samples.


Taxonomic binning
-----------------

In addition to taxonomic profiling, you can also predict the taxonomic identity
of your metagenomic samples by taxonomic binning. You can perform the taxonomic
binning with DNA-to-protein classifiers like `Kaiju <http://kaiju.binf.ku.dk/>`_
or with DNA-to-DNA classifiers like `Kraken2 <https://github.com/DerrickWood/kraken2/wiki>`_.

Using Kaiju
...........

First, let's do it through :code:`kaiju` command. To do so, we have
to choose which database we want Kaiju to download. In this case, we will only
focus on the viral communities of the metagenomic samples. Let's run the
:code:`kaiju` command like so:

.. code-block:: bash

    metabiome kaiju -i decontaminated_reads/ -o kaiju_out/ -x -k -d viruses

From this running, you will find two main output directories in the directory
:file:`kaiju_out/`: :file:`taxa_names/` and :file:`krona/`, which contain
the taxa classification of the assigned reads and their visualization through
krona figures, respectively.

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
    wget -P kraken2_db https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20201202.tar.gz
    tar -xvzf kraken2_db/k2_viral_20201202.tar.gz -C kraken2_db/

Now that we have a database, we can perform the taxonomic classification using
the following command:

.. code-block:: bash

    metabiome kraken2 -i decontaminated_reads/ -o kraken2_out/ -db kraken2_db/

Visualizing Kraken results
''''''''''''''''''''''''''

We have just performed the taxonomic classification of our reads with Kraken, so
let's visualize these results using `Krona <https://github.com/marbl/Krona/wiki>`_:

.. code-block:: bash

    metabiome krona -i kraken2_out/ -o krona_out/

And that's all! Inside the :file:`krona_out/` folder you will now find the Krona
graphs displaying the composition of your samples. Your result should be
similar to `this <_static/taxonomy.krona.html>`_.

Functional profiling
--------------------

The first time you use `HUMAnN <https://huttenhower.sph.harvard.edu/humann/>`_,
you must download two databases, ChocoPhlAn and a translated search database
(UniRef), see `HUMAnN documentation <https://github.com/biobakery/humann#5-download-the-databases>`_
for more info about this. Here we will download the demo version of ChocoPhlAn
database and the demo version of UniRef90 database by running the following
commands:

.. code-block:: bash

    # Activate environment containing HUMAnN
    conda activate metabiome-taxonomic-profiling

    # Create folder in which databases will be saved
    mkdir humann_db

    # Download databases
    humann_databases --download chocophlan DEMO humann_db/
    humann_databases --download uniref DEMO_diamond humann_db/

    # Deactivate environment
    conda deactivate

After downloading databases we are ready to profile our samples with HUMAnN:

.. code-block:: bash

    metabiome humann3 -i decontaminated_reads/ -o humann_results/


Extract 16S rDNA sequences
--------------------------
Now, lets suppose you want to perform additional analyses based on the 16S rDNA.
The :code:`bbduk` command can extract the 16S rDNA from your metagenomic samples through
`BBDuk <https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/>`_.
But first, you will need to download the 16S rDNA sequences from the database of
your choice. In this case, we will use our `custom 16S rDNA database of the phylum Firmicutes
<https://drive.google.com/file/d/1dOIgupiE-xpORIR-7jxaTMI63NXQBvdH/view?usp=sharing>`_.
Go ahead and run :code:`bbduk` command like so:


.. code-block:: bash

    metabiome bbduk -i decontaminated_reads/ -o bbduk_out/ \
        -D Firmicutes_rRNA_16S_silva.fa.gz -opts -Xmx2g

The output of :code:`bbduk` command is located in :file:`bbduk_out/`. This output is
very similar to the `Decontamination section <Decontamination_>`_ output.
However, in this context these files are the metagenomic reads that did
aligned to the Firmicutes 16S rDNA sequences.

*De-novo* Assembly
******************

Genome assembly
---------------

In this step you can use two different assemblers that receive the output from
:code:`bowtie2`: `metaSPAdes <https://cab.spbu.ru/software/spades/>`_ and
`MEGAHIT <https://github.com/voutcn/megahit>`_, in order to obtain contigs.
You can use just the assembler you like the most, or use both as we will do in
this tutorial. To perform the assembly, just run the following commands but keep
present that this may take several minutes so just sit tight!


Using MetaSPAdes
................

.. code-block:: bash

    # metaSPAdes
    metabiome metaspades -i decontaminated_reads/ -o metaspades_assembled_reads/


Using MEGAHIT
.............

.. code-block:: bash

    # MEGAHIT
    metabiome megahit -i decontaminated_reads/ -o megahit_assembled_reads/

.. note::

    By default, Metabiome doesn't perform co-assembly of multiple samples but
    instead it runs individual assemblies for each sample. If you want to
    perform co-assembly of many samples, see :ref:`How to perform co-assembly of
    samples <co-assembly>`.

These output genome draft assemblies are frequently used to perform genome quality assessment
and binning.

Quality assembly
----------------

In order to assess the quality of the assemblies performed in the previous step,
we are going to use `MetaQUAST <http://quast.sourceforge.net/metaquast>`_. The
minimal input for MetaQUAST is a folder with contigs in FASTA format, then
MetaQUAST will search and download reference sequences for you. However, in this
tutorial we will use the Metabiome's ``-opts`` flag (See :ref:`opts-flag`) in
order to give MetaQUAST a reference sequence to compare our contigs. As BeAn
58058 virus was one of the most abundant virus in our samples, we will use its
genome:

.. code-block:: bash

    # Create directory with reference sequence
    mkdir metaquast_ref_seq

    # Download reference genome
    wget -P metaquast_ref_seq ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/BeAn_58058_virus/latest_assembly_versions/GCF_001907825.1_ViralProj357638/GCF_001907825.1_ViralProj357638_genomic.fna.gz

    # Run MetaQUAST
    metabiome metaquast -i megahit_assembled_reads/ERR981212_sub_paired_bt2/ -o metaquast_out \
        -opts -r metaquast_ref_seqs/GCF_001907825.1_ViralProj357638_genomic.fna.gz

Genome binning
**************

The following step is to generate bins from the previous draft genomes or
contigs (wether from MetaSPAdes or MEGAHIT). To do so, we will use three
different binners: `MetaBAT2 <https://bitbucket.org/berkeleylab/metabat/>`_,
`MaxBin2 <https://sourceforge.net/projects/maxbin2/>`_ and `CONCOCT <https://concoct.readthedocs.io/en/latest/>`_.
Depending on the options you provide, these binners will need the contigs and
the reads that generated those contigs in order to run. In this case, we will
use both files located in the directory :file:`contigs_reads/`.

.. note:: Keep in mind that your contigs must have the same filename as
    their respective paired-end reads. Thus, your :file:`contigs_reads/`
    directory should look like this:

    .. code-block:: bash

        # Contig and their respective paired-end reads of the sample ERR981212
        ERR981212_sub_paired_bt2.fasta
        ERR981212_sub_paired_bt2_1.fq.gz
        ERR981212_sub_paired_bt2_2.fq.gz

        # Contig and their respective paired-end reads of the sample ERR981213
        ERR981213_sub_paired_bt2.fasta
        ERR981213_sub_paired_bt2_1.fq.gz
        ERR981213_sub_paired_bt2_2.fq.gz


Using MetaBAT2
--------------

Let's begin with MetaBAT2, which requires the contigs in gzip format in
order to run. Here is an example of how you should do it before running
:code:`metabat2` command :

.. code-block:: bash

    # Create input directory
    mkdir gzip_contigs

    # Copy contigs to the input directory
    cp contigs_reads/*.fasta gzip_contigs/

    # Compress the contigs in the required gzip format
    gzip gzip_contigs/*.fasta

    # Run MetaBAT2
    metabiome metabat2 -i gzip_contigs/ -o metabat2_out/ \
        -opts -m 1500 --maxP 50 --minS 30 --maxEdges 100 --minClsSize 1000

For example, MetaBAT2 will generate 23 bins from the assembly of
the sample ERR981212, which are located in :file:`metabat2_out/ERR981212_sub_paired_bt2/`.

.. code-block:: bash

    ERR981212_sub_paired_bt2.1.fa
    ERR981212_sub_paired_bt2.2.fa
    ERR981212_sub_paired_bt2.3.fa
    ERR981212_sub_paired_bt2.4.fa
    ......
    ERR981212_sub_paired_bt2.21.fa
    ERR981212_sub_paired_bt2.22.fa
    ERR981212_sub_paired_bt2.23.fa

Using MaxBin2
-------------

The next binner will be MaxBin2. Let's run the command :code:`maxbin2`
like so:

.. code-block:: bash

    metabiome maxbin2 -i contigs_reads/ -o maxbin2_out/ \
        -opts -min_contig_length 500 -prob_threshold 0.6

For example, MaxBin2 will generate just 1 bin and many too-short
bins from the sample ERR981212, which are located in
:file:`maxbin2_out/ERR981212_sub_paired_bt2/` and
:file:`maxbin2_out/ERR981212_sub_paired_bt2/ERR981212_sub_paired_bt2.tooshort`,
respectively.

Using CONCOCT
-------------

Last but not least, let's run :code:`concoct` command like so:

.. code-block:: bash

    metabiome concoct -i contigs_reads/ -o concoct_out/ -opts --no_original_data

For example, CONCOCT will generate 8 bins from the assembly
of the sample ERR981212, which are located in
:file:`concoct_out/fasta_bins/ERR981212_sub_paired_bt2/`:

.. code-block:: bash

    0.fa
    1.fa
    2.fa
    .....
    7.fa
    8.fa

.. note::

    In order to boost the binning process, you can also generate
    read-based coverage files that will help improve the bins,
    see :ref:`How to create read-based coverage files for genome
    binning <boost_binning>`.
