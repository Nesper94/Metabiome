.. _usermanual:

User manual
===========

.. contents::

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
