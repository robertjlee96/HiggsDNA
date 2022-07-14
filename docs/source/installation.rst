Installation
============

Standard Procedure
------------------

The standard way to get HiggsDNA up and running consists in the following steps.

After cloning the repo and accessing it with the usual::

        git clone https://gitlab.cern.ch/HiggsDNA-project/HiggsDNA.git
        cd HiggsDNA

one can create a conda environment with the main needed dependencies::

        conda env create -f environment.yml
        # if available, do it with mamba, it's much faster :)
        mamba env create -f environment.yml

        conda activate higgs-dna

To install the package in editable mode run::

        pip install -e .[dev]
