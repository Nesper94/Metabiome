.. _howto:

How-to...
=========

Run software natively
---------------------
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

.. _co-assembly:

...Perform co-assembly of samples
---------------------------------
