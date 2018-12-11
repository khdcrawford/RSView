==============================
Documentation for ``RSView``
==============================

The documentation is written in `reStructuredText` and can be built via `sphinx`.

Building the documentation
----------------------------
To build the documentation, you will need to have installed:

  * `sphinx`

This is included in the `rsview` conda environment.

Once the environment is properly setup and acitvated, and the `rsview`package  has been installed by running `python setup.py install` within the ``RSView`` directory, you can build the documentation with::

    make html

The HTML documentation will then be in ``./_build/html/``.
