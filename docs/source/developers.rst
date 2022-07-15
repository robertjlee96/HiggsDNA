Instructions for Developers
===========================

Please lint, format and run tests before sending a PR:

.. code-block:: bash

   flake8 higgs_dna
   black higgs_dna
   pytest

--------------------
Update Documentation
--------------------

When the package is modified, it is highly recommended to also update the documentation accordingly. The files that make the docs that you are reading are located inside ``docs/source``. You can modify them or add new ones.

When you are satisfied with the result, do the following:

.. code-block:: bash

   cd docs
   sphinx-apidoc -o source/modules ../higgs_dna
   sphinx-build source build/html

The ``sphinx-apidoc`` command will build the documentation from the package's docstrings, so every change in the package itself will be picked up.
At this point you can see (locally) how the updated docs look like by simply opening the just built html section using your favourite browser, e.g.:

.. code-block:: bach

   firefox build/html/index.html
