Main Concepts
=============

.. _def-cltool:

Command Line Tool
-----------------

If you want to run an analysis with processors and taggers that have already been developed, the suggested way is to use the command line tool ``run_analysis.py``.

.. _def-processor:

Processor
---------
Processors are items defined within Coffea where the analysis workflow is described. While a general overview is available in the `Coffea documentation <https://coffeateam.github.io/coffea/concepts.html#coffea-processor>`_, here we will focus on the aspects that are important for HiggsDNA.

Since in Higgs to diphoton analysis there are some operations that are common to every analysis workflow, we wrote a base processor which can be used in many basics analysis. If more complex operations are needed, one can still write a processor that inherits from the base class, but rewriting the function `process`. The operations that one can find within HggBaseprocessor.process are the following:

- application of filters and triggers
- Chained Quantile Regression to correct shower shapes and isolation variables
- photon IdMVA
- diphoton IdMVA
- photon preselection
- event tagging
- application of systematic uncertainties
