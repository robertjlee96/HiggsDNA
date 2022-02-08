=============
Main Concepts
=============


.. _def-cltool:

-----------------
Command Line Tool
-----------------
If you want to run an analysis with processors and taggers that have already been developed, the suggested way is to use the command line tool ``run_analysis.py``.

Following are the main parts of the analysis that can be configured through the command line:

* datasets
  path to a JSON file in the form ``{"dataset_name": [list_of_files]}`` (like the one dumped by dasgoclient)
* workflow
  the coffea processor you want to use to process you data, can be found in the modules located inside the subpackage ``higgs_dna.workflows``
* metaconditions
  the name (without ``.json`` extension) of one of the JSON files inherited from FLASHgg and located inside ``higgs_dna.metaconditions``
* taggers
  the set of taggers you want to use, can be found in the modules located inside the subpackage ``higgs_dna/workflows/taggers``

Here is how these information `can` be passed to through command line:

.. code-block:: bash

        run_analysis.py --samples DY-MC.json --workflow tagandprobe --meta Era2017_legacy_xgb_v1 --tagger-set DummyTagger1

However, when the analysis increases in complexity, the parameters specified can increase a lot in complexity and the command can get quite messy (see Combine as a perfect example of this). In these cases, the above mentioned parameters can be specified in a JSON file and passed to the command line with the flag ``--json-analysis``

.. code-block:: bash

        run_analysis.py --json-analysis simple_analysis.json

where ``simple_analysis.json`` looks like this:

.. code-block:: json

        {
            "samplejson": "/work/gallim/devel/HiggsDNA/tmp/DY-data-test.json",
            "workflow": "tagandprobe",
            "metaconditions": "Era2017_legacy_xgb_v1",
            "taggers": [
                "DummyTagger1"
            ]
        }


.. _def-processor:

---------
Processor
---------
Processors are items defined within Coffea where the analysis workflow is described. While a general overview is available in the `Coffea documentation <https://coffeateam.github.io/coffea/concepts.html#coffea-processor>`_, here we will focus on the aspects that are important for HiggsDNA.

Since in Higgs to diphoton analysis there are some operations that are common to every analysis workflow, we wrote a base processor which can be used in many basics analysis. If more complex operations are needed, one can still write a processor that inherits from the base class and redefines the function ``process``. The operations that one can find within ``HggBaseprocessor.process`` are the following:

* application of filters and triggers
* Chained Quantile Regression to correct shower shapes and isolation variables
* photon IdMVA
* diphoton IdMVA
* photon preselection
* event tagging
* application of systematic uncertainties
