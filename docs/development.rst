Development
===========

Run tests
---------

Use the ``nencarta`` conda environment from the repository root:

.. code-block:: powershell

   $env:NUMBA_CACHE_DIR = "$PWD\ignore\_pytest_cache\numba"
   $env:NENCARTA_CACHE_DIR = "$PWD\ignore\_pytest_cache\nencarta"
   python -m pytest -q -p no:cacheprovider

The test suite includes:

* configuration validation tests
* workspace path tests
* ``process_watershed`` entrypoint tests with the pipeline mocked
* local path checks for the South Platte calls used in ``ignore/new.py`` and
  ``ignore/og.py``
* a SHA-256 manifest check for cached output folders in ``ignore/``

Refresh ignored-output baselines
--------------------------------

The ignored-output baseline does not copy large raster or vector files. It stores
relative filenames, file sizes, and SHA-256 hashes in
``tests/fixtures/ignore_output_manifest.json``.

Regenerate the manifest only after you intentionally accept new outputs:

.. code-block:: powershell

   python tests\write_ignore_output_manifest.py

The manifest currently covers these local output folders:

* ``ignore/forecast_test_original``
* ``ignore/forecast_test_new``
* ``ignore/clean_dem_original``
* ``ignore/clean_dem_new``
* ``ignore/no_bathy_original``
* ``ignore/no_bathy_new``

When refactoring, run the tests before and after the code change. The manifest
test must remain identical unless the intended change is to alter generated
products.

Build the documentation locally
-------------------------------

Install the documentation dependencies:

.. code-block:: bash

   pip install -r docs/requirements.txt

Then build the HTML output from the repository root:

.. code-block:: bash

   sphinx-build -b html docs docs/_build/html

Read the Docs configuration
---------------------------

This repository is configured for Read the Docs with:

* ``.readthedocs.yaml`` at the repository root
* ``docs/conf.py`` for the Sphinx project configuration
* ``docs/requirements.txt`` for documentation-only Python dependencies

When the repository is connected to a Read the Docs project, the platform will build the documentation from the ``docs/`` directory automatically.
