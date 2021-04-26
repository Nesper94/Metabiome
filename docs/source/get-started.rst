.. _getstarted:

Get started
===========

.. contents::

Requirements
------------

Metabiome runs in Linux systems and requires the `Conda
<https://conda.io/miniconda.html>`_ environment manager to manage the software
packages it uses.

Installing Miniconda
^^^^^^^^^^^^^^^^^^^^
The simplest way to get Conda is by installing `Miniconda
<https://docs.conda.io/en/latest/miniconda.html>`_.
In order to install Miniconda please follow the `Miniconda instructions
<https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_
and be sure to execute :code:`conda init` at the end of the installation
process.

Installing Metabiome
--------------------

In order to install, clone Metabiome into your machine and execute
:file:`install.sh` by running the following commands:

.. code-block:: bash

    git clone -b master https://github.com/Nesper94/metabiome.git
    cd metabiome/
    bash install.sh

You can also download Metabiome directly from the `Github page
<https://github.com/Nesper94/metabiome>`_.

The installation script will create Conda environments with all the software
you can use with Metabiome.

For a more complete reference see the :ref:`user manual <usermanual>`.