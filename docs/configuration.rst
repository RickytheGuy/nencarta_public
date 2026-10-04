Configuration Reference
=======================

JSON Inputs
--------------

Below is a description of all the JSON inputs that NenCarta accepts when you run the
``json`` subcommand. Each watershed definition is a JSON object with the required
keys and optional processing parameters described above. You can include multiple
watershed objects in the ``watersheds`` array to run them in batch mode.

.. _json-age_of_forecast_days:

* ``age_of_forecast_days`` (Integer, optional): This is how old previous forecasts
  will be allowed to be based upon the current day. For example, specifying 7 days
  will delete all forecast flood inundation maps that are older than 7 days. Default
  7.

.. _json-bathy_args:

* ``bathy_args`` (Dictionary, optional): A dictionary of arguments that will be
  passed to ARC and the flood mapper when estimating bathymetry in each cross-section.
  See the ARC and the respective flood mapper documentation for more information on
  these arguments. See the `ARC documentation <https://github.com/MikeFHS/automated-rating-curve/wiki/ARC-Input-File-Arguments>`_ 
  for more information on the arguments that can be passed to ARC. See the 
  `Curve2Flood documentation <https://github.com/MikeFHS/curve2flood/wiki/Curve2Flood-Input-File-Arguments>`_
  for more information on the arguments that can be passed to Curve2Flood.

.. _json-bathy_use_banks:

* ``bathy_use_banks`` (Bool, optional): Setting this to True will allow the system
  to use estimated bank elevations to estimate bathymetry using either the
  streamflow specified in ``specified_bathyflow_field`` or the optional
  drainage-area power-law depth relationship. Setting this value to False will
  allow the system to use the estimated water surface elevation to estimate
  bathymetry using that same bathymetry target. Default False.

.. _json-clean_dem:

* ``clean_dem`` (Bool, optional): Setting this value to True will allow the system to
  use the DEM Cleaner. Setting this value to False will bypass use of the DEM Cleaner.
  Default False.

.. _json-create_reach_average_curve_file:

* ``create_reach_average_curve_file`` (Bool, optional): Setting this value to true
  instructs NenCarta to direct ARC to output the curve file by taking velocity, depth,
  and top-width estimates for all stream cells on an individual stream reach and
  creating the curve parameters for each stream cell. Setting this value to False
  instructs NenCarta to direct ARC to output the curve file using velocity, depth, and
  top-width estimates for each stream cell. Default False.

.. _json-dem_dir:

* ``dem_dir`` (String): A full filepath to the directory containing one or more DEMs
  that you will be using as input in NenCarta. If 
  :ref:`move_stream_network_to_thalweg <json-move_stream_network_to_thalweg>`
  is ``true``, and a DEM is stored in geographic coordinates
  (latitude/longitude), NenCarta first creates a projected GeoTIFF copy in the
  watershed ``FlowDirection`` folder using ``WGS 84 / NSIDC EASE-Grid 2.0
  Global`` (EPSG:6933). That projected copy becomes the DEM used for all later
  processing for that tile.

.. _json-dem_filter:

* ``dem_filter`` (String, optional): A glob string with which files in the
  ``dem_dir`` must match to be included in the run. By default, "*", or all files.

.. _json-disable_bathymetry:

* ``disable_bathymetry`` (Bool, optional): Setting this equal to true will disable
  the estimation of bathymetry within the hydraulic calculations and the creation of a
  composite topobathymetric surface.

.. _json-estimate_consequences:

* ``estimate_consequences`` (Bool, optional): Setting this equal to true will utilize
  Go-Consequences and the National Structure Inventory (NSI) to perform consequence
  assessment for the area within each flood inundation map. This functionality is
  currently only available for the coterminous United States. Default False.

.. _json-find_banks_based_on_landcover:

* ``find_banks_based_on_landcover`` (Bool, optional): Setting this value to True will
  direct NenCarta to first try finding the banks using the land cover. If the stream
  cell is in a water pixel, NenCarta will direct ARC to search for the end of the
  water surrounding the stream cell and this will location will be designated as the
  banks. Setting this value to False will direct NenCarta and ARC to find the banks by
  assuming the channel is flat in the DEM and allowing ARC to find the banks when the
  flat water surface ends in the cross-section. Default True.

.. _json-flood_waterlc_and_strm_cells:

* ``flood_waterlc_and_strm_cells`` (Bool, optional): Argument for the flood mapper,
  whether to force in the flood map all cells that are wet in the land cover and
  stream raster as wet in the flood map. Default False.

.. _json-floodmap_args:

* ``floodmap_args`` (Dictionary, optional): A dictionary of arguments that will be
  passed to the flood mapper when creating flood inundation maps. See the respective
  documentation for more information on these arguments. See the 
  `ARC input-file documentation <https://github.com/MikeFHS/automated-rating-curve/wiki/ARC-Input-File-Arguments>`_ 
  for more information on the arguments that can be passed to ARC. See the 
  `Curve2Flood documentation <https://github.com/MikeFHS/curve2flood/wiki/Curve2Flood-Input-File-Arguments>`_
  for more information on the arguments that can be passed to Curve2Flood.

.. _json-floodmap_identifier:

* ``floodmap_identifier`` (String, optional): A string that will be appended to the
  flood map filenames to help identify the flood maps created by this NenCarta
  simulation. Default is an empty string.

.. _json-floodmap_mode:

* ``floodmap_mode`` (String, optional): Either ``forecast``, ``user``, or
  ``return_period``. Forecast mode runs GEOGLOWS or NWM forecast flows. User
  mode uses ``user_flow_files``. Return-period mode uses ``return_periods``.
  Default ``forecast``.

.. _json-flowline:

* ``flowline`` (String): A full filepath to the flowline shapefile that you'll be
  using to run NenCarta.

.. _json-forensic_forecast_date:

* ``forensic_forecast_date`` (String, optional): If you want to use a past forecast,
  you can specify the date here in YYYYMMDD format (e.g., "20250807" is August 7,
  2025) in UTC. The archive of GEOGLOWS forecasts goes back to July 1, 2024. We are
  unsure of the current limitations for the National Water Model archive of forecasts.

.. _json-forensic_forecast_hour:

* ``forensic_forecast_hour`` (String, optional): If you want to use a past National
  Water Model forecast, along with the forensic_forecast_date, you must also specify a
  forensic forecast hour. This is the hour the forecast was produced, in UTC. For the
  "NWM_short_range" forecast, the forensic forecast hour must be between 0 and 23,
  expressed as a string (e.g., "00"). For the "NWM_medium_range" forecast, the
  forensic forecast hour must be "00", "06", "12", or "18". For the
  "NWM_long_range" forecast, the forensic forecast hour must be "00".

.. _json-geoglows_vpu:

* ``geoglows_vpu`` (Integer, optional): Only used for CONUS Vector Processing Units
  (VPUs) in GEOGLOWS. When the NenCarta user wants to use the GEOGLOWS map-tables to
  identify where flooding is forecast in CONUS, the GEOGLOWS VPU can be specified in
  this option. NenCarta only requires that the streamline vector for the VPU be
  available in this use case. NenCarta will use the location of flowlines forecasts to
  meet or exceed the 2-year discharge to both download USGS 3DEP DEM data and simulate
  flood inundation.

.. _json-lake_filter_json:

* ``lake_filter_json`` (String, optional): The path to the GEOGLOWS json that
  describes which stream reaches are within a lake. Optional input. This is currently
  only functional for GEOGLOWS data.

.. _json-land_watervalue:

* ``land_watervalue`` (Int, optional): The value in the land cover raster that
  represents water. Default 80.

.. _json-make_ap_database:

* ``make_ap_database`` (Bool, optional): Whether ARC will make an Area-Perimeter
  file. Default True.

.. _json-make_curvefile:

* ``make_curvefile`` (Bool, optional): Whether ARC will make a Curve file. Default
  True.

.. _json-make_depth_maps:

* ``make_depth_maps`` (Bool, optional): Whether or not to make flood depth maps.
  Default True.

.. _json-make_fist_inputs:

* ``make_fist_inputs`` (Bool, optional): Whether or not to make FIST inputs for
  flood. Default True.

.. _json-make_velocity_maps:

* ``make_velocity_maps`` (Bool, optional): Whether or not to make flood velocity
  maps. Default True.

.. _json-make_wse_maps:

* ``make_wse_maps`` (Bool, optional): Whether or not to make flood water surface
  elevation maps. Default True.

.. _json-mannings_text_file:

* ``mannings_text_file`` (String, optional): The full filepath to the Manning's n
  text file to be used in the flood mapping process. If not provided, NenCarta will
  use its built-in default Manning's n values.

.. _json-mapper:

* ``mapper`` (String, optional): Here you're specifying if you're running
  "FloodSpreader", "Curve2Flood-Kernel Weighted", or "Curve2Flood-FLDPLNpy", or
  "Curve2Flood-Multi-Point Interpolation" when performing bathymetry estimation and
  flood inundation mapping. Defaults to FloodSpreader

.. _json-move_stream_network_to_thalweg:

* ``move_stream_network_to_thalweg`` (Bool, optional): If True, NenCarta moves
  each reach of the source flowlines onto the D8 flow path it follows in the DEM,
  keeping its ID, so that every reach is one unbroken flow path (see
  :ref:`stream-network-movement`). Default is False. This is required if using
  FLDPLNpy as the ``mapper`` option.

.. _json-name:

* ``name`` (String): The name of the watershed you're modeling.

.. _json-new_strm_threshold_km2:

* ``new_strm_threshold_km2`` (Float, optional): The smallest drainage area, in
  square km, of a DEM flow path that source reaches are moved onto. Where a source
  reach runs through cells draining less than this (the upper end of a
  headwater, typically), that part of it is left out. Rivers flowing into the
  domain from outside it count from where they enter. Required when
  ``move_stream_network_to_thalweg`` is True.

.. _json-num_workers:

* ``num_workers`` (Integer, optional): Worker count used when ``parallel`` is True
  and NenCarta creates its own process executor. If omitted, NenCarta uses the
  default executor behavior.

.. _json-nwm_api_key:

* ``nwm_api_key`` (String, optional): Required when streamflow_source is set to "NWM"
  or any "NWM_*" option. This is the NWM API key passed as the 'x-api-key' header for
  NWM requests. You must apply for an API key using `these instructions 
  <https://docs.ciroh.org/docuhub-staging/docs/products/data-management/bigquery-api/>`_.

.. _json-output_dir:

* ``output_dir`` (String): The full filepath to the directory where your output will
  be saved.

.. _json-folder_paths:

* ``folder_paths`` (Dictionary, optional): Where to put each of the output
  folders, as a path template per folder. The folders are ``DEM``,
  ``ARC_InputFiles``, ``FloodMap``, ``Bathymetry``, ``DEM_Updated``, ``STRM``,
  ``LAND``, ``FLOW``, ``VDT``, ``ESA_LC``, ``FIST``, ``Consequences``,
  ``FlowDirection`` and ``FLDPLN``; a ``default`` entry applies to every folder
  not listed.
  Templates can use ``{output_dir}``, ``{name}`` (the watershed name), ``{dem}``
  (the DEM's file name without its extension, or the bounding box name) and
  ``{folder}`` (the folder's name), and relative paths are relative to
  ``output_dir``. For example, ``{"FloodMap": "/results/FloodMap/{dem}",
  "VDT": "/results/VDT/{dem}"}`` sends the flood maps and VDTs of every DEM to
  one tree and leaves everything else where it was, and ``{"default":
  "{folder}/{name}"}`` groups every folder by type instead of by watershed. Keep
  ``{name}`` or ``{dem}`` in a template when several watersheds or DEMs share
  ``output_dir``: some files, such as ``fldpln_library.parquet``, are not named per
  DEM. File names do not change (see ``short_file_names`` and ``file_names``).
  Default None, which puts every folder in ``{output_dir}/{name}/{folder}``.

.. _json-short_file_names:

* ``short_file_names`` (Bool, optional): Name files by what they are alone,
  dropping the ``<NWM|GEOGLOWS>_<DEM>_`` or ``<DEM>_`` prefix, e.g.
  ``VDT/VDT_Database_Bathy.parquet`` instead of
  ``VDT/GEOGLOWS_<DEM>_buffered_VDT_Database_Bathy.parquet``. Since the DEM is no
  longer in the names, every DEM needs its own folders: NenCarta raises an error
  if two DEMs would share one, as they do by default when ``dem_dir`` holds
  several DEMs, in which case put ``{dem}`` in ``folder_paths``. Buffered flood
  maps are cropped in place rather than renamed. Default False.

.. _json-file_names:

* ``file_names`` (Dictionary, optional): New names for kinds of files, keyed by
  the part of the name that says what the file is, e.g. ``{"VDT_Database_Bathy":
  "vdt"}``; the prefix and extension are kept, so with ``short_file_names`` that
  file is ``VDT/vdt.parquet``. Files that come one per return period, forecast or
  flow file keep that part after the new name, e.g. ``{"ARC_Flood": "flood"}``
  names the return period flood maps ``flood_rp2.tif`` ... ``flood_rp100.tif``
  (with ``short_file_names``), and an empty name names them by that part alone,
  e.g. ``{"ARC_Flood": ""}`` names them ``rp2.tif`` ... ``rp100.tif``; a file of
  that kind without such a part keeps its usual name. Files that differ only by
  extension, like the ``matched`` stream vector and raster, are renamed
  together. The keys are
  ``dem``, ``fixed``, ``filled``, ``fixed_Clean``, ``Clean``, ``flowdir``,
  ``StrmShp``, ``STRM_Raster``,
  ``STRM_Raster_Clean``, ``matched``, ``lakes``, ``stream_info``,
  ``fldpln_library``, ``LAND_Raster``, ``AR_Manning_n_MED``, ``Reanalysis``,
  ``2yr_flow_initial``, ``rp`` (the return period flow files, ``rp<N>.csv``),
  ``forecast``, ``Flow_COMID_Q``, ``ARC_Input_Bathy``, ``ARC_Input_InitialFlood``,
  ``ARC_Input_FloodForecast``, ``ARC_Input`` (the flood mapper configs),
  ``VDT_Database_Initial``, ``VDT_Database_Bathy``, ``AP_Database_Bathy``,
  ``CurveFile``, ``CurveFile_Initial``, ``CurveFile_Bathy``, ``XS``,
  ``Representative_XS``, ``water_mask``, ``ARC_Bathy``, ``FS_Bathy``,
  ``ARC_Flood``, ``ARC_Flood_Initial``, ``ARC_Flood_Bathy``, ``ARC_Depth``,
  ``ARC_FloodDepth``, ``ARC_FloodWSE``, ``ARC_FloodVEL``, ``Seed`` and ``FIST``
  (the FIST GeoJSON files). Names cannot contain ``/`` or ``\``, and NenCarta
  raises an error if two files would end up with the same path. Default None.

.. _json-bbox:

* ``bbox`` (List of Numbers, optional): Bounding box as
  ``[minx, miny, maxx, maxy]``. Used with ``source_dems`` when NenCarta builds a
  DEM for the domain.

.. _json-dem:

* ``dem`` (String, optional): Path to a single DEM file. ``dem_dir`` takes
  precedence when both ``dem`` and ``dem_dir`` are supplied.

.. _json-source_dems:

* ``source_dems`` (List of Strings, optional): Source DEM files used to create a
  domain DEM when ``bbox`` is supplied.

.. _json-source_flowlines:

* ``source_flowlines`` (List of Strings, optional): Source flowline files used
  when NenCarta prepares the flowline network instead of using one prepared
  ``flowline`` file.

.. _json-buffer:

* ``buffer`` (Bool, optional): Whether to buffer the requested ``bbox`` when
  creating a DEM from ``source_dems``. Default False.

.. _json-buffer_distance:

* ``buffer_distance`` (Float, optional): Buffer distance applied when ``buffer``
  is True. Default 0.1.

.. _json-land_cover_cache:

* ``land_cover_cache`` (List of Strings, optional): Local land-cover raster
  cache files. Default empty list.

.. _json-lakes:

* ``lakes`` (String, optional): Path to a lake polygon file. With
  ``burn_streams``, lakes are left out of the water mask the burned DEM is
  smoothed within and the bathymetry is bounded by, and when
  ``move_stream_network_to_thalweg`` is True they are cut out of the moved
  network. With the ``Curve2Flood-FLDPLNpy`` mapper no stream runs inside a lake:
  a reach flowing into one ends at its shore, and a reach crossing one keeps its
  longest stretch outside it. With the other mappers a reach touching a lake is
  dropped unless it carries a river through it, with reaches above and below it
  that touch no lake.

.. _json-reanalysis_file:

* ``reanalysis_file`` (String, optional): Existing streamflow reanalysis file to
  use instead of generating or downloading one.

.. _json-return_period_file:

* ``return_period_file`` (String, optional): Return period flows to build the
  reanalysis file from, instead of the GEOGLOWS retrospective return periods. A
  local path or a URL (``s3://``, ``https://``) to a CSV, Parquet, Zarr (``.zarr``)
  or NetCDF (``.nc``) file. Zarr and NetCDF files need ``river_id`` and
  ``return_period`` dimensions; CSV and Parquet files need a ``river_id`` column
  and one ``rp<N>`` column per return period (``rp2``, ``rp5``, ...). The river
  ID may instead be named after ``stream_id_field``. A NetCDF file on S3 or HTTP is
  downloaded once per process, so Zarr is faster for large remote datasets.
  GEOGLOWS only. Default None.

.. _json-return_period_variables:

* ``return_period_variables`` (List of Strings, optional): Variables of
  ``return_period_file`` that hold return period flows. When there are several,
  the highest is used for each river and return period. Default None, which uses
  every variable with exactly the ``river_id`` and ``return_period`` dimensions:
  ``gumbel``, ``gumbel_hourly`` and ``gumbel_daily`` for GEOGLOWS.

.. _json-include_fdc:

* ``include_fdc`` (Bool, optional): Whether to add flow duration curve
  exceedance flows (``p_exceed_0``, ``p_exceed_5``, ... ``p_exceed_100`` and
  ``p_exceed_1``) to the reanalysis file. GEOGLOWS only. Default True.

.. _json-fdc_file:

* ``fdc_file`` (String, optional): Flow duration curves to take the exceedance
  flows from, instead of the GEOGLOWS retrospective flow duration curves. Read
  like ``return_period_file``, with a ``p_exceed`` dimension (exceedance
  probability in percent) for Zarr and NetCDF, or ``p_exceed_<P>`` columns for
  CSV and Parquet. It needs all of the exceedances listed under ``include_fdc``.
  Unlike the GEOGLOWS default, which falls back to computing the curves from the
  daily retrospective flows when it cannot be read, an ``fdc_file`` that cannot be
  read is an error. Default None.

.. _json-fdc_variable:

* ``fdc_variable`` (String, optional): Variable of the flow duration curve
  dataset to use. Default None, which uses ``hourly_annual`` if the dataset has it
  (as GEOGLOWS does), or else the dataset's only variable with exactly the
  ``river_id`` and ``p_exceed`` dimensions.

.. _json-reanalysis_storage_options:

* ``reanalysis_storage_options`` (Dictionary, optional): fsspec storage options
  for a remote ``return_period_file`` or ``fdc_file``, e.g.
  ``{"profile": "my-aws-profile"}`` for a private S3 bucket. Default None, which
  reads S3 anonymously.

.. _json-raise_errors_if_river_ids_missing:

* ``raise_errors_if_river_ids_missing`` (Bool, optional): What to do when some
  river IDs in the domain are not in the return period or flow duration curve
  dataset, e.g. for a domain that crosses the edge of a regional dataset. When
  True, raise an error; when False, log a warning and leave those rivers out of
  the reanalysis file, so it only has the rivers found in every dataset used.
  Default True.

.. _json-return_periods:

* ``return_periods`` (List of Integers, optional): Return periods to map when
  ``floodmap_mode`` is ``return_period``. Valid values are 2, 5, 10, 25, 50,
  and 100.

.. _json-parallel:

* ``parallel`` (Bool, optional): Whether to run workspace tasks in parallel.
  Default False.

.. _json-profile:

* ``profile`` (Bool, optional): Whether to print task profiling output for
  serial runs. Default True.

.. _json-compression:

* ``compression`` (String, optional): Raster compression option passed to output
  writers. Default ``LZW``.

.. _json-exclude:

* ``exclude`` (List of Integers, optional): Stream IDs to exclude from the run.
  Default empty list.

.. _json-use_parquet:

* ``use_parquet`` (Bool, optional): Whether to use parquet intermediate files
  where supported. Default False.

.. _json-use_yaml:

* ``use_yaml`` (Bool, optional): Whether to write mapper configuration files as
  YAML instead of text. Default False.

.. _json-use_vrt:

* ``use_vrt`` (Bool, optional): Whether to create assigned DEMs as VRT files
  where supported. Default False.

.. _json-streams_as_parquet:

* ``streams_as_parquet`` (Bool, optional): Whether to write prepared stream
  vectors as parquet files. Default False.

.. _json-burn_streams:

* ``burn_streams`` (Bool, optional): Whether to burn the source streams into the
  DEM, writing ``DEM_Updated/{DEM}_fixed.tif``. Each reach is lowered along its
  whole line, in steps of 0.5 m that never rise downstream, upstream reaches
  before the reaches they flow into, and where reaches meet the larger one keeps
  the shared cell, so the burned channels drain continuously and filling the DEM
  does not flatten them. Default False.

.. _json-project_to_utm:

* ``project_to_utm`` (Bool, optional): Whether to project outputs to a local UTM
  coordinate system where supported. Default False.

.. _json-raise_errors_if_nothing_in_domain:

* ``raise_errors_if_nothing_in_domain`` (Bool, optional): Whether empty-domain
  checks should raise errors instead of logging and continuing. Default True.

.. _json-use_power_laws_for_bathymetry:

* ``use_power_laws_for_bathymetry`` (Bool, optional): Whether to use drainage
  area power laws for bathymetry parameters when supported. Default False.

.. _json-use_dem_derived_channel_mask:

* ``use_dem_derived_channel_mask`` (Bool, optional): Whether stream burning uses
  a channel mask derived from half-meter DEM elevations instead of the existing
  bathymetry water mask. Used by the stream burning path in ``move_streams``.
  Default False.

.. _json-area_m2_field:

* ``area_m2_field`` (String, optional): Drainage area field in square meters.
  Default ``DSContArea``.

.. _json-area_km2_field:

* ``area_km2_field`` (String, optional): Drainage area field in square
  kilometers. Default ``DSContArea_km2``.

.. _json-coefficient_depth:

* ``coefficient_depth`` (Float, optional): Coefficient for power-law bathymetry
  depth estimation when ``use_power_laws_for_bathymetry`` is enabled. Default
  0.27.

.. _json-exponent_depth:

* ``exponent_depth`` (Float, optional): Exponent for power-law bathymetry depth
  estimation when ``use_power_laws_for_bathymetry`` is enabled. Default 0.21.

.. _json-coefficient_width:

* ``coefficient_width`` (Float, optional): Coefficient for power-law bathymetry
  width estimation when ``use_power_laws_for_bathymetry`` is enabled. Default
  2.44.

.. _json-exponent_width:

* ``exponent_width`` (Float, optional): Exponent for power-law bathymetry width
  estimation when ``use_power_laws_for_bathymetry`` is enabled. Default 0.34.

.. _json-make_vdt:

* ``make_vdt`` (Bool, optional): Whether to write a VDT database. Default True.

.. _json-make_cross_section_file:

* ``make_cross_section_file`` (Bool, optional): Whether to write cross-section
  output from ARC. Default False.

.. _json-make_representative_cross_section_file:

* ``make_representative_cross_section_file`` (Bool, optional): Whether to write
  ARC's representative cross sections, which describe each reach by its cross
  sections averaged stage by stage, every 0.10 m, for routing in hydrologic models.
  ARC builds them in its bathymetry run, alongside the rating curves and from the
  same cross sections, so they cover the same reaches and carved channels as the
  VDT database, which they leave unchanged. Default False.

.. _json-fldpln_dh:

* ``fldpln_dh`` (Float, optional): FLDPLN depth interval. Default 0.5.

.. _json-fldpln_min_depth:

* ``fldpln_min_depth`` (Float, optional): Minimum FLDPLN depth. Default 0.1.

.. _json-fldpln_max_depth:

* ``fldpln_max_depth`` (Float, optional): Maximum FLDPLN depth. Default 25.0.

.. _json-fldpln_max_drop_below_source:

* ``fldpln_max_drop_below_source`` (Float, optional): Leave out of the FLDPLN
  library the floodplain cells more than a segment's maximum depth plus this
  many meters below the stream pixel that floods them. FLDPLN charges nothing
  for a downhill step, so every stream pixel reaches the whole slope beneath
  it, and on alluvial fans and bajadas those cells are most of
  the library. The flood mapper discards them anyway, keeping a cell only
  while it is no more than the stream pixel's depth plus
  ``FLDPLN_Max_Drop_Below_Source`` (in ``floodmap_args``) below it, and that
  depth comes from the same VDT as the segment's maximum depth, smoothed along
  the stream. Set this a few meters more than ``FLDPLN_Max_Drop_Below_Source``
  and the maps do not change: on a desert tile, 4 m more cut the library from
  213 million rows to 10 million, and peak memory from 18.7 GB to 2.3 GB, with
  every return period map identical, while equal values changed under 0.1% of
  the wet cells. Cells left out do not count toward
  ``fldpln_max_cells_per_segment``. Default None, which keeps them all.

.. _json-fldpln_max_cells_per_segment:

* ``fldpln_max_cells_per_segment`` (Integer, optional): The most floodplain
  cells (FLDPLN library rows) one stream segment may contribute. On truly flat
  ground, such as coastal plains, a segment's floodplain spreads level with
  the stream rather than below it, so ``fldpln_max_drop_below_source`` cannot
  shrink it, and a few segments can reach millions of cells each. A segment
  that would exceed this stops at the depth where it fills up instead, exactly
  as if its maximum depth were that depth, so its maximum depth is found from
  the terrain rather than from the VDT; the cells it floods at no depth at all
  are always kept. This bounds memory but is lossy: floods deeper than that
  depth are mapped only as far as it. On a coastal plain tile pruned with
  ``fldpln_max_drop_below_source``, 1,000,000 cut peak memory from 6.3 GB to
  4.6 GB, leaving the 2 and 5 year maps unchanged but losing 0.9% of the 50
  year and 4.3% of the 100 year wet area; 500,000 cut it to 3.4 GB and lost
  5.2% and 8.4%. Needs ``fldpln_solver`` ``exact``. Default None, which lets
  every segment flood to its full maximum depth.

.. _json-fldpln_max_wse_rise:

* ``fldpln_max_wse_rise`` (Float, optional): Maximum WSE rise used by FLDPLN.
  Default 0.01.

.. _json-fldpln_keep_spilling:

* ``fldpln_keep_spilling`` (Bool, optional): Whether FLDPLN keeps spilling
  cells. Default False.

.. _json-fldpln_parallel:

* ``fldpln_parallel`` (Bool, optional): Whether FLDPLN runs its internal parallel
  mode. Default False.

Runtime cache settings
----------------------

NenCarta can use local caches to avoid repeated metadata and compiled-function
work.

* ``NENCARTA_CACHE_DIR``: Optional directory for the LMDB raster metadata cache.
  If it is not set, NenCarta uses ``Path.home() / ".cache"``. If LMDB is not
  installed, the cache cannot be opened, or the directory is locked, NenCarta
  disables this cache and continues.
* ``NUMBA_CACHE_DIR``: Optional directory for Numba caches used by ARC and other
  compiled functions. This is useful on systems where the default user cache is
  locked or unavailable.

.. _json-overwrite_floodmaps:

* ``overwrite_floodmaps`` (Bool, optional): Whether to replace existing flood
  indunation, depth, WSE, etc. maps if they exist. Defaults to True.

.. _json-overwrite:

* ``overwrite`` (Bool, optional): Whether to rebuild existing intermediate
  products, including DEM-specific stream geometry, stream rasters, flow files,
  bathymetry inputs, and model setup files. Setting this value to False reuses
  existing intermediate products when they are present. Default False.

.. _json-q_baseflow_threshold:

* ``q_baseflow_threshold`` (Float, optional): Setting this value (in cubic meters per
  second) equal to a float value will filter out all streams in your domain that have
  a baseflow (from ``specified_bathyflow_field``) that are less than the specified
  value. Default is None.

.. _json-slope_low_percentile:

* ``slope_low_percentile`` (Float, optional): Lower percentile used when ARC
  derives robust stream slope values. Default 5.

.. _json-slope_high_percentile:

* ``slope_high_percentile`` (Float, optional): Upper percentile used when ARC
  derives robust stream slope values. Default 95.

.. _json-quiet:

* ``quiet`` (Bool, optional): Setting this value to True will suppress ARC and
  Curve2Flood output. Default False.

.. _json-remove_old_forecast_files:

* ``remove_old_forecast_files`` (Bool, optional): If True, check existing flood maps
  and see if they have forecast data older than ``age_of_forecast_days``. If so,
  remove. Default True.

.. _json-specified_bathyflow_field:

* ``specified_bathyflow_field`` (String, optional): The field in the GEOGLOWS
  downloaded reanalysis data that will be provided to ARC to estimate bathmetry in
  each cross-section. For "GEOGLOWS" it must be one of"p_exceed_0", "p_exceed_5",
  "p_exceed_10", "p_exceed_15", "p_exceed_20", "p_exceed_25", "p_exceed_30",
  "p_exceed_35", "p_exceed_40", "p_exceed_45", "p_exceed_50", "p_exceed_65",
  "p_exceed_70", "p_exceed_75", "p_exceed_80", "p_exceed_85", "p_exceed_90",
  "p_exceed_95", "p_exceed_100", "rp2", "rp5", "rp10", "rp25", "rp50",
  "rp100","p_exceed_0_premium", or "rp100_premium". For "NWM" it must be one of
  "rp2", "rp5", "rp10", "rp25", "rp50", "rp100", or "rp100_premium". Default is
  "p_exceed_50". If the five drainage-area bathymetry parameters are supplied,
  NenCarta omits ``Flow_File_BF`` from the ARC input file and uses the
  drainage-area relationship for ARC bathymetry instead. In that case,
  ``specified_bathyflow_field`` is only still relevant if you also use
  ``q_baseflow_threshold`` to filter streams before the ARC run.

.. _json-specify_depths_for_bathy_mask:

* ``specify_depths_for_bathy_mask`` (List of Numbers, optional): If
  ``use_specified_depth_for_bathy_mask`` is True, the user must specify at least one
  depth value for FloodSpreaderPy to use when creating a water mask for the DEM
  cleaner and the bathymetry estimation processes. If ``clean_dem`` is False, this
  argument requires only one float value and that value will be used to create the
  bathymetry estimation water mask. If ``clean_dem`` is True, you will need to
  specify two float values, the first value is the depth used to create the water mask
  for the DEM cleaner and the second value is used to create the bathymetry estimation
  water mask. The json processes accepts these inputs as a Python list. The command
  line interface accepts these as one or two separate values (i.e.,
  --specify_depths_for_bathy_mask 1.0 2.0).

.. _json-specified_highflow_field:

* ``specified_highflow_field`` (String, optional): The field in the GEOGLOWS
  downloaded reanalysis data that will be provided to ARC as the highest flow used to
  estimate water surface elevation. For "GEOGLOWS" it must be one of "p_exceed_0",
  "p_exceed_5", "p_exceed_10", "p_exceed_15", "p_exceed_20", "p_exceed_25",
  "p_exceed_30", "p_exceed_35", "p_exceed_40", "p_exceed_45", "p_exceed_50",
  "p_exceed_65", "p_exceed_70", "p_exceed_75", "p_exceed_80", "p_exceed_85",
  "p_exceed_90", "p_exceed_95", "p_exceed_100", "rp2", "rp5", "rp10", "rp25",
  "rp50", "rp100","p_exceed_0_premium", or "rp100_premium". For "NWM" it must be
  one of "rp2", "rp5", "rp10", "rp25", "rp50", "rp100", or "rp100_premium".
  Default is "rp100_premium".

.. _json-streamflow_source:

* ``streamflow_source`` (String, optional): Setting this equal to "GEOGLOWS" will
  force NenCarta to use GEOGLOWS retrospective and forecast streamflow data. The
  deafult is GEOGLOWS. Setting this to "NWM_short_range" will force NenCarta to use
  the National Water Model retrospective and short-range forecast streamflow data.
  Setting this to "NWM_medium_range" will force NenCarta to use the National Water
  Model retrospective and medium-range forecast streamflow data. Setting this to
  "NWM_long_range" will force NenCarta to use the National Water Model retrospective
  and long-range forecast streamflow data. Default "GEOGLOWS".

.. _json-StrmOrder_Field:

* ``StrmOrder_Field`` (String, optional): The field in the flowline GIS data that
  specifies the stream order of the streams in your model domain. This input is
  required if you plan to use StrmOrder_Lower or StrmOrder_Upper to limit which
  streams will be used for flood inundation mapping by NenCarta.

.. _json-StrmOrder_Lower:

* ``StrmOrder_Lower`` (Integer, optional): The lowest value of stream order that you
  plan to use in your NenCarta simulation.

.. _json-StrmOrder_Upper:

* ``StrmOrder_Upper`` (Integer, optional): The highest value of stream order that you
  plan to use in your NenCarta simulation.

.. _json-use_specified_depth_for_bathy_mask:

* ``use_specified_depth_for_bathy_mask`` (Bool, optional): Setting this value to
  False will direct NenCarta to create a water mask that is a compilation of the
  stream network raster and the land use designated as water. This will be used to
  clean the input DEM and burn bathymetry into the DEM. Setting this value to True
  will direct NenCarta to create a water mask for DEM cleaning that fills the stream
  cell with a specified depth of water. The depth of water will be based upon the
  value(s) specified in the ``specify_depths_for_bathy_mask`` argument. Default True.

.. _json-use_warning_flags_to_download_dem:

* ``use_warning_flags_to_download_dem`` (Bool, optional): Whether to use warning
  flag locations to identify DEM downloads for the selected GEOGLOWS VPU.
  Default False.

.. _json-user_flow_files:

* ``user_flow_files`` (String or List of Strings, optional): If ``floodmap_mode`` is
  "user", than use this file or list of files to create flood maps instead of looking
  at the forecast.

.. _json-vdt_file_extension:

* ``vdt_file_extension`` (String, optional): The file extension of the VDT files to
  be created. Default "txt".

Required watershed keys
-----------------------

Each watershed definition requires these arguments in the JSON file. They are defined in 
`this section <https://nencarta.readthedocs.io/en/latest/configuration.html#json-inputs>`_.

* :ref:`name <json-name>`
* :ref:`flowline <json-flowline>`
* :ref:`dem_dir <json-dem_dir>`
* :ref:`output_dir <json-output_dir>`

Common processing options
-------------------------

These options control the main flood-mapping workflow. They are defined in 
`this section <https://nencarta.readthedocs.io/en/latest/configuration.html#json-inputs>`_.

* :ref:`overwrite <json-overwrite>`
* :ref:`clean_dem <json-clean_dem>`
* :ref:`use_warning_flags_to_download_dem <json-use_warning_flags_to_download_dem>`

Forecast and flow options
-------------------------

These options select the streamflow source and forecast behavior. They are defined in 
`this section <https://nencarta.readthedocs.io/en/latest/configuration.html#json-inputs>`_.

* :ref:`streamflow_source <json-streamflow_source>`: Forecast provider and product selection.
* :ref:`forensic_forecast_date <json-forensic_forecast_date>`: Use a prior forecast by date in ``YYYYMMDD`` format.
* :ref:`forensic_forecast_hour <json-forensic_forecast_hour>`: Required for archived NWM forecasts.
* :ref:`age_of_forecast_days <json-age_of_forecast_days>`: Remove outdated forecast outputs based on age.
* :ref:`remove_old_forecast_files <json-remove_old_forecast_files>`: Enable cleanup of stale forecast products.
* :ref:`user_flow_files <json-user_flow_files>`: Provide user-supplied flow files when :ref:`floodmap_mode <json-floodmap_mode>` is ``user``.

ARC options
-----------
These options control the ARC workflow that estimates bathymetry and creates curve files. They are defined in 
`this section <https://nencarta.readthedocs.io/en/latest/configuration.html#json-inputs>`_.

* :ref:`create_reach_average_curve_file <json-create_reach_average_curve_file>`
* :ref:`make_curvefile <json-make_curvefile>`
* :ref:`make_ap_database <json-make_ap_database>`
* :ref:`make_representative_cross_section_file <json-make_representative_cross_section_file>`
* :ref:`bathy_args <json-bathy_args>`
* :ref:`vdt_file_extension <json-vdt_file_extension>`

ARC outputs
~~~~~~~~~~~

For each DEM tile, NenCarta writes bathymetry-related outputs under
``output_dir`` / ``name``: 

* ``VDT/<NWM|GEOGLOWS>_<DEM>_VDT_Database_Bathy.<vdt_file_extension>``: Bathymetry
  VDT database.
* ``VDT/<NWM|GEOGLOWS>_<DEM>_VDT_FS_Bathy.csv``: Bathymetry VDT test file.
* ``VDT/<NWM|GEOGLOWS>_<DEM>_CurveFile_Bathy.csv``: Bathymetry curve file, written
  when ``make_curvefile`` is ``true``.
* ``VDT/<NWM|GEOGLOWS>_<DEM>_AP_Database_Bathy.txt``: Area-perimeter database,
  written when ``make_ap_database`` is ``true``.
* ``VDT/<NWM|GEOGLOWS>_<DEM>_Representative_XS.<csv|parquet>``: Representative
  cross sections, one row per reach and 0.10 m stage, written by the ARC
  bathymetry run when ``make_representative_cross_section_file`` is ``true``. It
  is Parquet when ``vdt_file_extension`` is ``parquet`` and CSV otherwise.


Bathymetry options
------------------

These options control bathymetry estimation generation. They are defined in 
`this section <https://nencarta.readthedocs.io/en/latest/configuration.html#json-inputs>`_.

* :ref:`disable_bathymetry <json-disable_bathymetry>`
* :ref:`bathy_use_banks <json-bathy_use_banks>`
* :ref:`drainage_area_field <json-drainage_area_field>`
* :ref:`coefficient_depth <json-coefficient_depth>`
* :ref:`exponent_depth <json-exponent_depth>`
* :ref:`coefficient_width <json-coefficient_width>`
* :ref:`exponent_width <json-exponent_width>`
* :ref:`use_specified_depth_for_bathy_mask <json-use_specified_depth_for_bathy_mask>`
* :ref:`specify_depths_for_bathy_mask <json-specify_depths_for_bathy_mask>`
* :ref:`find_banks_based_on_landcover <json-find_banks_based_on_landcover>`

Bathymetry outputs
~~~~~~~~~~~~~~~~~~

For each DEM tile, NenCarta writes bathymetry-related outputs under
``output_dir`` / ``name``. ARC and the flood mapper (e.g., Curve2Flood) 
create the main bathymetry products:

* ``Bathymetry/<NWM|GEOGLOWS>_<DEM>_ARC_Bathy.tif``: ARC bathymetry raster written
  from the ARC bathymetry run.
* ``Bathymetry/<NWM|GEOGLOWS>_<DEM>_FS_Bathy.tif``: FloodSpreader or Curve2Flood
  bathymetry raster used later as the DEM input for flood mapping when bathymetry is
  enabled.
* ``FloodMap/<NWM|GEOGLOWS>_<DEM>_ARC_Flood[_floodmap_identifier]_Bathy.tif``:
  Bathymetry flood raster generated during the mapper bathymetry stage.
* ``FloodMap/<NWM|GEOGLOWS>_<DEM>_ARC_Flood[_floodmap_identifier]_Bathy.shp``:
  Bathymetry flood geometry written alongside the bathymetry flood raster.

If :ref:`disable_bathymetry <json-disable_bathymetry>` is ``true``, the code still writes the hydraulic support
files in ``VDT/`` needed for later flood mapping, but it returns before creating the
``ARC_Bathy.tif``, ``FS_Bathy.tif``, and bathymetry flood-map outputs.


Flood-map options
-----------------
These options control the flood-map generation process. They are defined in 
`this section <https://nencarta.readthedocs.io/en/latest/configuration.html#json-inputs>`_.

* :ref:`mapper <json-mapper>`
* :ref:`make_depth_maps <json-make_depth_maps>`
* :ref:`make_velocity_maps <json-make_velocity_maps>`
* :ref:`make_wse_maps <json-make_wse_maps>`
* :ref:`floodmap_identifier <json-floodmap_identifier>`
* :ref:`floodmap_args <json-floodmap_args>`

Flood-map outputs
~~~~~~~~~~~~~~~~~

Flood-map outputs are written under ``output_dir`` / ``name`` / ``FloodMap``. The
exact filenames vary depending on whether you are using forecast flows or
:ref:`user_flow_files <json-user_flow_files>`:

* Flood extent raster: NenCarta always writes a flood extent raster through
  ``OutFLD``. For forecast mode this becomes a forecast-named file such as
  ``<NWM|GEOGLOWS>_<DEM>_ARC_Flood[_floodmap_identifier]_Forecast_<date>[ _<hour>].tif``.
  For user-flow mode the flow-file stem is appended to the flood-map filename.
* Flood geometry: NenCarta writes a companion geometry file through ``OutSHP`` when
  :ref:`floodmap_args <json-floodmap_args>` does not set ``Make_Output_GPKG`` to
  ``False``. The path uses the
  same base name as the flood extent raster and changes the suffix to ``.shp``.
* Flood depth raster: ``<streamflow_source>_<DEM>_ARC_FloodDepth[_floodmap_identifier]...``
  is written when :ref:`make_depth_maps <json-make_depth_maps>` is ``true``.
* Flood water-surface-elevation raster:
  ``<streamflow_source>_<DEM>_ARC_FloodWSE[_floodmap_identifier]...`` is written when
  :ref:`make_wse_maps <json-make_wse_maps>` is ``true``.
* Flood velocity raster:
  ``<streamflow_source>_<DEM>_ARC_FloodVEL[_floodmap_identifier]...`` is written when
  :ref:`make_velocity_maps <json-make_velocity_maps>` is ``true``.

The flood mapper uses ``FS_Bathy.tif`` as the DEM input when bathymetry is enabled.
If bathymetry is disabled, it falls back to the cleaned DEM when
:ref:`clean_dem <json-clean_dem>` is ``true`` or the original DEM otherwise.

In forecast mode, NenCarta produces one set of flood outputs per DEM tile for the
selected forecast. In user-flow mode, NenCarta produces one set of flood outputs per
input flow file.

Consequence estimation
----------------------

Set :ref:`estimate_consequences <json-estimate_consequences>` to ``true`` to invoke the Go-Consequences workflow after flood-map generation.
See `this section <https://nencarta.readthedocs.io/en/latest/configuration.html#json-inputs>`_ for more information on 
what ``estimate_consequences`` does.

Consequence estimation outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When :ref:`estimate_consequences <json-estimate_consequences>` is ``true``, ``nencarta/main.py`` iterates over the
depth rasters returned by the flood-mapping stage and creates one consequences run per
depth raster. The outputs are:

* ``FloodMap/<depth_file_stem>_WGS84.tif``: A WGS84 copy of the depth raster created
  with ``gdal.Warp`` before the Docker run. If reprojection fails, the original depth
  raster is used instead.
* ``Consequences/<depth_file_stem>_WGS84_consequences.json`` or
  ``Consequences/<depth_file_stem>_consequences.json``: The Go-Consequences config
  JSON written by ``Create_Go_Consequence_GeoJSON``.
* ``Consequences/<depth_file_stem>_WGS84_consequences.gpkg`` or
  ``Consequences/<depth_file_stem>_consequences.gpkg``: The GeoPackage results file
  written by the Go-Consequences Docker container.

The consequences JSON points Go-Consequences to the flood-depth raster in
``/data/FloodMap`` and writes the results GeoPackage to ``/data/Consequences`` inside
the container, which corresponds to the watershed output directory on the host.

If :ref:`estimate_consequences <json-estimate_consequences>` is ``true`` and
:ref:`make_depth_maps <json-make_depth_maps>` is ``false``, the code forces
:ref:`make_depth_maps <json-make_depth_maps>` back to ``true`` before processing so
that a depth
raster exists for the consequences workflow.

FIST options
------------

NenCarta can create inputs for the Flood Inundation Surface Topology (FIST) model as a forecast. To do this, set
``make_fist_inputs`` to ``true`` to generate FIST-ready inputs as part of the workflow.
NenCarta writes those files to the ``FIST`` subdirectory in the
``output_dir``/ ``name`` directory.

NenCarta relies on ARC to generate the FIST inputs. For the underlying ARC workflow, see
the `ARC FIST documentation <https://automated-rating-curve.readthedocs.io/en/latest/making_inputs_for_fist/>`_.


FIST outputs
~~~~~~~~~~~~

The generated FIST inputs are created for each forecast scenario, including minimum,
median, and maximum streamflow forecasts. They include stream-cell locations, water
surface elevations, and SEED values (0 = not a seed, 1 = a seed). SEED values designate 
the furthest upstream points for headwater streams and are used by FIST to define 
flow paths and inundation extents.

The FIST subdirectory contains the following file types:

* ``*_{forecast_date}_min.geojson``: GeoJSON point features for stream-cell locations,
  with water-surface elevation and SEED values for the minimum streamflow forecast.

* ``*_{forecast_date}_med.geojson``: GeoJSON point features for stream-cell locations,
  with water-surface elevation and SEED values for the median streamflow forecast.

* ``*_{forecast_date}_max.geojson``: GeoJSON point features for stream-cell locations,
  with water-surface elevation and SEED values for the maximum streamflow forecast.

* ``*_Seed.shp``: A shapefile containing the SEED locations that designate the furthest
  upstream points for headwater streams.

.. _stream-network-movement:

Stream network movement options
-------------------------------

* :ref:`move_stream_network_to_thalweg <json-move_stream_network_to_thalweg>`
* :ref:`new_strm_threshold_km2 <json-new_strm_threshold_km2>`
* :ref:`burn_streams <json-burn_streams>`
* :ref:`lakes <json-lakes>`

FLDPLN, ARC and Curve2Flood read the stream network off the DEM's D8 flow
directions: FLDPLN walks each reach down the D8 pointers, so a reach that wanders
off its flow path stops being one reach. Source flowlines (TDX-Hydro, GEOGLOWS,
NHD, ...) were drawn from other DEMs and drift on and off this DEM's flow paths.
When ``move_stream_network_to_thalweg`` is enabled, NenCarta moves each source
reach onto the D8 path it follows, keeping its ID so that the flows keyed by it
still apply. The work is done in ``nencarta/tasks/move_streams.py`` and
``nencarta/tasks/stream_conflation.py``.

1. With ``burn_streams``, the source streams are burned into the DEM first, so
   that its flow paths follow them.
2. `WhiteboxTools <https://github.com/jblindsay/whitebox-tools>`_ fills the
   DEM's depressions and derives its D8 flow directions.
3. Every source reach is sampled about once a cell, from its upstream end down,
   and matched onto the D8 network as a hidden Markov model: the candidates at
   each sample are the nearby cells draining at least ``new_strm_threshold_km2``,
   a candidate can only be followed by one downstream of it, and the match keeps
   the stretch of D8 path that best follows the line. A reach the DEM doesn't
   follow at all is left out rather than forced somewhere wrong.
4. The matched reaches are laid onto the D8 network in flow order, each running
   from its matched start to where the next reach starts. Where two reaches' flow
   paths meet somewhere the source network has no confluence, the source network
   decides which carries on; a reach whose match ends short of the reach it flows
   into carries on down its D8 path to meet it.

Because each reach is matched on its own, a reach the DEM disagrees with costs
only that reach, never the basin above it. Reaches are split exactly where the
source network splits them, including the breaks that are not confluences. An
endorheic source reach ends where its line ends, even though the filled DEM
routes its flow on over the basin's spill point. Each moved reach's downstream ID
is the reach its last cell drains into on the D8 network, or -1, and every reach
is one unbroken D8 path.

If ``overwrite`` is ``false`` and the moved stream network products already
exist, NenCarta reuses them instead of rebuilding them.


Stream network movement outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``DEM_Updated/{DEM}_fixed.tif``: the DEM with the source streams burned in,
  when ``burn_streams`` is enabled.

* ``DEM_Updated/{DEM}_filled.tif``: the depression-filled DEM the flow
  directions are derived from.

* ``FlowDirection/{DEM}_flowdir.tif``: the D8 flow direction raster, in
  WhiteboxTools' pointer codes.

* ``STRM/{DEM}_matched.gpkg`` (``.parquet`` with ``streams_as_parquet``): the
  moved stream network, used by every later step. Each row is a source reach
  that was placed: its ID, the ID of the reach it drains into (-1 for none), a
  ``topological_order`` that ranks every reach before the one it drains into,
  and the source reach's other attributes. Its line runs through the centres of
  its D8 cells, from upstream to downstream.

* ``STRM/{DEM}_matched.tif``: the moved stream network's cells, each holding its
  reach's ID.

* ``STRM/{DEM}_lakes.tif``: the rasterized ``lakes``, when given.

* ``FLDPLN/{DEM}_stream_info.parquet`` (``.csv`` without ``use_parquet``): each moved
  reach's D8 walk for the FLDPLN mapper.


GUI options
-----------

The GUI is a runner and monitor for ``process_watershed``. It validates the
form inputs, previews the watershed dictionary that will be passed to the Python
API, starts the run, and streams log output while the simulation is active. It
does not provide map visualization.

The GUI inputs are organized into tabs that correspond to the categories of JSON
inputs described above.

Required Inputs
~~~~~~~~~~~~~~~

* ``Watershed Name`` -> ``name``
* ``Flowline File`` -> ``flowline``
* ``Source Flowline Files`` -> ``source_flowlines``
* ``DEM Directory`` -> ``dem_dir``
* ``DEM File`` -> ``dem``
* ``Source DEM Files`` -> ``source_dems``
* ``Bounding Box`` -> ``bbox``
* ``Output Directory`` -> ``output_dir``
* ``Folder Path Templates`` -> ``folder_paths``
* ``Short File Names`` -> ``short_file_names``
* ``File Names`` -> ``file_names``

Key Workflow Switches
~~~~~~~~~~~~~~~~~~~~~

* ``Clean DEM (Requires Initial Flood Map Step)`` -> ``clean_dem``
* ``Estimate Consequences (Run Go-Consequences)`` -> ``estimate_consequences``
* ``Mapper Method`` -> ``mapper``
* ``Streamflow Source`` -> ``streamflow_source``
* ``NWM API Key`` -> ``nwm_api_key``
* ``Parallel`` -> ``parallel``
* ``Worker Count`` -> ``num_workers``

For ``Streamflow Source``, the GUI labels map to the JSON values as follows:

* ``GEOGLOWS`` -> ``GEOGLOWS``
* ``NWM Short Range`` -> ``NWM_short_range``
* ``NWM Medium Range`` -> ``NWM_medium_range``
* ``NWM Long Range`` -> ``NWM_long_range``

Advanced Parameters
~~~~~~~~~~~~~~~~~~~

* ``Disable Bathymetry`` -> ``disable_bathymetry``
* ``Bathy Use Banks`` -> ``bathy_use_banks``
* ``Flood LC and Stream Cells in Flood Map`` -> ``flood_waterlc_and_strm_cells``
* ``Use Specified Depth for Bathy Mask`` -> ``use_specified_depth_for_bathy_mask``
* ``Find Banks Based on Land Cover (Default=True)`` -> ``find_banks_based_on_landcover``
* ``Overwrite Existing Products`` -> ``overwrite``
* ``Create Reach Average Curve File`` -> ``create_reach_average_curve_file``
* ``Use Warning Flags to Download DEM`` -> ``use_warning_flags_to_download_dem``
* ``Land Water Value`` -> ``land_watervalue``
* ``Age of Forecast Days`` -> ``age_of_forecast_days``
* ``Specific flood depths (in meters) for bathy mask`` -> ``specify_depths_for_bathy_mask``
* ``GEOGLOWS VPU ID`` -> ``geoglows_vpu``
* ``Forensic Forecast Date`` -> ``forensic_forecast_date``
* ``Forensic Forecast Hour`` -> ``forensic_forecast_hour``
* ``Bathy Flow Field`` -> ``specified_bathyflow_field``
* ``High Flow Field`` -> ``specified_highflow_field``
* ``Use DEM-derived Channel Mask`` -> ``use_dem_derived_channel_mask``
* ``Use Power Laws for Bathymetry`` -> ``use_power_laws_for_bathymetry``
* ``Depth Coefficient`` -> ``coefficient_depth``
* ``Depth Exponent`` -> ``exponent_depth``
* ``Width Coefficient`` -> ``coefficient_width``
* ``Width Exponent`` -> ``exponent_width``
* ``Move Stream Network to Thalweg`` -> ``move_stream_network_to_thalweg``
* ``Stream Threshold for New Stream Network`` -> ``new_strm_threshold_km2``
* ``Slope Low Percentile`` -> ``slope_low_percentile``
* ``Slope High Percentile`` -> ``slope_high_percentile``
* ``Stream Order Field`` -> ``StrmOrder_Field``
* ``Stream Order Lower`` -> ``StrmOrder_Lower``
* ``Stream Order Upper`` -> ``StrmOrder_Upper``
* ``Baseflow Threshold`` -> ``q_baseflow_threshold``
* ``Return Period File`` -> ``return_period_file``
* ``Return Period Variables`` -> ``return_period_variables``
* ``Include Flow Duration Curve Exceedances`` -> ``include_fdc``
* ``Flow Duration Curve File`` -> ``fdc_file``
* ``Flow Duration Curve Variable`` -> ``fdc_variable``
* ``Error On River IDs Missing From Flow Datasets`` -> ``raise_errors_if_river_ids_missing``
* ``Lake Filter JSON`` -> ``lake_filter_json``
* ``Overwrite Forecast Floodmaps`` -> ``overwrite_floodmaps``
* ``Remove Old Forecast Files`` -> ``remove_old_forecast_files``
* ``Make FIST Inputs`` -> ``make_fist_inputs``
* ``DEM Filter`` -> ``dem_filter``
* ``Floodmap Mode`` -> ``floodmap_mode``
* ``User Flow Files`` -> ``user_flow_files``
* ``Make Curve File`` -> ``make_curvefile``
* ``Make Area-Perimeter Database`` -> ``make_ap_database``
* ``Make Representative Cross-Section File`` -> ``make_representative_cross_section_file``
* ``Make Depth Maps`` -> ``make_depth_maps``
* ``Make Velocity Maps`` -> ``make_velocity_maps``
* ``Make WSE Maps`` -> ``make_wse_maps``
* ``VDT File Extension`` -> ``vdt_file_extension``
* ``Manning's n Text File`` -> ``mannings_text_file``
* ``Pre-processsing/Bathymetry Arguments`` -> ``bathy_args``
* ``Floodmap Arguments`` -> ``floodmap_args``

GUI-specific conversions
~~~~~~~~~~~~~~~~~~~~~~~~

* The GUI stores ``Watershed Name`` as ``watershed_name`` internally and then
  sends it to ``process_watershed`` as ``name``.
* ``Specific flood depths (in meters) for bathy mask`` is entered as comma-separated
  text in the GUI and written as a JSON list for
  ``specify_depths_for_bathy_mask``.
* Multi-file controls are sent as JSON lists.
