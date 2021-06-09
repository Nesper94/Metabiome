.. _howto:

How-to...
=========

...Run software natively
------------------------
If you need to run software in a way Metabiome does not support, you can run it
natively by simply activating the environment where it lives and then using it
normally. In the following example taken from the Bowtie 2 manual, we need to
use this software to create an index for the Lambda phage reference genome, and
not for decontamination of reads as implemented in Metabiome:

.. code-block:: bash

    # Activate environment
    conda activate metabiome-preprocessing
    # Run Bowtie 2 natively
    bowtie2-build $BT2_HOME/example/reference/lambda_virus.fa lambda_virus

You can find a complete list of Metabiome modules and software in :ref:`modules`.

.. _co-assembly:

...Perform co-assembly of samples
---------------------------------

If you need to perform an assembly using the reads from all or many samples at
once, then you can concatenate the files together and run the assemblers as if
it where just one paired-end sample. For example, assuming you have files named
according to the ``_1.``/``_2.`` naming convention:

.. code-block:: bash

    # Create directory to store new files
    mkdir co-assembly

    # Concatenate all _1. files
    cat *_1.fq.gz > co-assembly/all_1.fq.gz

    # Concatenate all _2. files
    cat *_2.fq.gz > co-assembly/all_2.fq.gz

    # Run assembly
    metabiome metaspades -i co-assembly -o assembly_results

.. _boost_binning:

...Create read-based coverage files for genome binning
------------------------------------------------------

If you want to create read-based coverage files to boost the
binning of your contigs, you can do it easily
through our pipeline Metabiome. Metabiome has three different
software (CONCOCT, MaxBin2 and MetaBAT2) for genome binning. So,
let's just take a look on how to create these read-based coverage
files and bin our contigs alongside:

For CONCOCT and MaxBin2 is pretty straightforward because the commands
:code:`concoct` and :code:`maxbin2` create the read-based coverage files
automatically like so:

.. code-block:: bash

    # Run CONCOCT
    metabiome concoct -i contigs_reads/ -o concoct_out/

    # Run MaxBin2
    metabiome maxbin2 -i contigs_reads/ -o maxbin2_out/

Whereas for MetaBAT2, if you want to perform the binning with read-based
coverage files, you must generate them previously. To do so you can use
our command :code:`coverm` like so:

.. code-block:: bash

    # Run CoverM
    metabiome coverm  -i contigs_reads/ -o read_coverage/

Now, let's use this read coverage files to run :code:`metabat2` command.
But keep in mind that MetaBAT2 requires the contigs in gzip format in
order to run:

.. code-block:: bash

    # Run MetaBAT2
    metabiome metabat2 -i gzip_contigs/ -co read_coverage/ -o metabat2_out/


Moreover, if you already have your read-based coverage files, you can also
provide it to any of the binners and automatically skip the process of
creating these read-based coverage files.

.. _scaffolds2bin:

...create scaffolds-to-bin tsv files for DAS Tool
-------------------------------------------------

Besides of contigs in fasta format, DAS Tool requires a tab separated
scaffolds-to-bin file from each metagenomic binner. To create these files,
you must have the bins from each metagenomic binner in distinct directories.
Now, we will use the bins from MaxBin2's output (:file:`maxbin_bins/`).

.. code-block:: bash

    # Move to directory containing bins
    cd ~/maxbin_bins/

    # Activate environment
    conda activate metabiome-das_tool

    # Create scaffolds-to-bin tsv files
    Fasta_to_Scaffolds2Bin.sh -e fasta > ERR981212_scaffolds2bin_maxbin2.tsv

    # Deactivate environment
    conda deactivate metabiome-das_tool

In this case, we have succesfully created scaffolds-to-bin tsv files
of MaxBin2's ouput from sample ERR981212. You can now do it with
other metagenomic binners. Keep in mind that, for example, MetaBAT2
generate bins with :code:`fa` extension. So, you must change the
:code:`-e` flag from :code:`fasta` to :code:`fa`.
