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
