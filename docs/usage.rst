Usage
=====

NenCarta supports three primary execution modes:

* GUI mode
* JSON-driven runs
* Direct CLI arguments

You can also call the Python entrypoint directly with
``nencarta.process_watershed`` when you already have a watershed dictionary.

GUI mode
--------

Launch the built-in graphical interface with:

.. code-block:: bash

   flood-mapping gui

The GUI is a runner and monitor. It validates a watershed configuration,
previews the dictionary that will be passed to ``process_watershed``, starts the
run, and streams log output while the simulation is active. It does not include
map visualization.

JSON mode
---------

Use the ``json`` subcommand when you want to run one or more watersheds from a structured configuration file.

Example:

.. code-block:: json

   {
     "watersheds": [
       {
         "name": "yellowstone_example",
         "flowline": "C:/path/to/flowline.shp",
         "dem_dir": "C:/path/to/dem_dir",
         "output_dir": "C:/path/to/output",
         "overwrite": true,
         "mapper": "FloodSpreader",
         "streamflow_source": "NWM_short_range",
         "nwm_api_key": "YOUR_NWM_API_KEY"
       }
     ]
   }

Run the file in serial mode:

.. code-block:: bash

   flood-mapping json "/path/to/your.json" --serial

Run the file in parallel mode:

.. code-block:: bash

   flood-mapping json "/path/to/your.json" --parallel --num_workers 8

CLI mode
--------

Use the ``cli`` subcommand to run a single watershed directly from the terminal:

.. code-block:: bash

   flood-mapping cli ExampleWatershed "C:\path\to\flowline.shp" "C:\path\to\dem_dir" "C:\path\to\output" --overwrite --mapper FloodSpreader --streamflow_source NWM_short_range --nwm_api_key "YOUR_NWM_API_KEY"

Python entrypoint
-----------------

``process_watershed`` is the public Python entrypoint used by the CLI, JSON
runner, GUI, and local scripts in ``ignore/``. It accepts one watershed
dictionary. The function validates the dictionary with ``NencartaConfig``,
creates one ``Workspace`` per matching DEM, and sends those workspaces through
the pipeline.

Example:

.. code-block:: python

   from nencarta import process_watershed

   process_watershed(
       {
           "name": "south_platte_example",
           "flowline": r"C:\path\to\streams.parquet",
           "dem_dir": r"C:\path\to\DEMs",
           "output_dir": r"C:\path\to\outputs",
           "dem_filter": "fabdem.tif",
           "mapper": "Curve2Flood-Kernel Weighted",
           "streamflow_source": "GEOGLOWS",
           "overwrite": True,
           "use_specified_depth_for_bathy_mask": True,
           "specify_depths_for_bathy_mask": [0.1],
       }
   )

DEM selection
~~~~~~~~~~~~~

NenCarta uses the following priority to decide which DEMs become workspaces:

* ``dem_dir`` plus ``dem_filter``: one workspace per matching DEM.
* ``dem``: one workspace for the supplied file.
* ``source_dems`` plus ``bbox``: one synthetic workspace that builds the DEM for
  the domain.
* No DEM source: the run logs a warning and returns without starting the
  pipeline.

If you want ARC bathymetry estimated from drainage-area power laws instead of a
baseflow field, provide all five optional parameters together:
``--drainage_area_field``, ``--coefficient_depth``, ``--exponent_depth``,
``--coefficient_width``, and ``--exponent_width``.

Forecast sources
----------------

NenCarta can use:

* ``GEOGLOWS``
* ``NWM_short_range``
* ``NWM_medium_range``
* ``NWM_long_range``

If you select an NWM source, you must supply ``nwm_api_key``.
