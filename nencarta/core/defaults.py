"""
Every configuration key nencarta accepts, with its default and a short description, grouped by topic.

The command line options are made from ``CONFIG_OPTIONS``, and the tests check that the GUI and the docs cover
every key, so a new key is added here first. ``DEFAULT_CONFIG`` holds just the defaults.
"""
from dataclasses import dataclass

from nencarta.core.enumerations import FloodMapMode, Mapper, StreamflowSource


VALID_RETURN_PERIODS = (2, 5, 10, 25, 50, 100)

# The flow columns of the reanalysis file nencarta builds
FLOW_FIELD_OPTIONS = (
    [f"p_exceed_{i}" for i in [0] + list(range(5, 101, 5))]
    + ["p_exceed_0_premium"]
    + ["rp2", "rp5", "rp10", "rp25", "rp50", "rp100", "rp100_premium"]
)

DEFAULT_BATHY_ARGS = {
    "VDT_Database_NumIterations": 30,
    "Make_Output_GPKG": True,
    "TW_MultFact": 1.5,
    "TopWidthPlausibleLimit": 2000,
    "Bathy_Trap_H": 0.2,
    "Degree_Manip": 6.1,
    "Degree_Interval": 1.5,
    "Low_Spot_Range": 2,
    "Gen_Dir_Dist": 10,
    "Gen_Slope_Dist": 10,
    "Stream_Slope_Method": "local_average_corrected",
}

DEFAULT_FLOODMAP_ARGS = {
    "Make_Output_GPKG": True,
    "TW_MultFact": 1.5,
    "TopWidthPlausibleLimit": 6000,
}


@dataclass(frozen=True)
class ConfigOption:
    """
    One configuration key. ``kind`` is the type of its value: text, path, int, float, bool, dict, or a list of
    text, paths, ints or floats (str_list, path_list, int_list, float_list). ``choices`` are its usual values.
    """
    default: object
    kind: str
    help: str
    choices: tuple = ()


CONFIG_SECTIONS = {
    "Required": {
        "name": ConfigOption(None, "text", "Name of the watershed. Its outputs go in output_dir/name."),
        "output_dir": ConfigOption(None, "path", "Folder to write the outputs in."),
    },
    "DEMs": {
        "dem": ConfigOption(None, "path", "A DEM file to map."),
        "dem_dir": ConfigOption(None, "path", "A folder of DEMs. Each DEM matching dem_filter is mapped on its own. "
                                              "Used instead of dem when both are given."),
        "dem_filter": ConfigOption("*", "text", "File pattern that picks the DEMs in dem_dir, e.g. *.tif."),
        "source_dems": ConfigOption([], "path_list", "DEMs to cut the domain's DEM from, with bbox or buffer."),
        "bbox": ConfigOption(None, "float_list", "Area to map, [minx, miny, maxx, maxy] in degrees. Its DEM is cut "
                                                 "from source_dems (and dem, if given)."),
        "buffer": ConfigOption(False, "bool", "Map a larger area, taken from source_dems, then crop the flood maps "
                                              "back to the DEM. Rivers then flow in from outside the DEM."),
        "buffer_distance": ConfigOption(0.1, "float", "How far buffer grows the DEM, in degrees."),
        "use_vrt": ConfigOption(False, "bool", "Write the domain's DEM as a VRT instead of a GeoTIFF."),
        "project_to_utm": ConfigOption(False, "bool", "Project the domain's DEM, when it is cut from source_dems, to "
                                                      "its UTM zone."),
        "use_warning_flags_to_download_dem": ConfigOption(False, "bool", "Download USGS DEMs into dem_dir where the "
                                                                         "GEOGLOWS forecast of forensic_forecast_date "
                                                                         "floods rivers in geoglows_vpu (CONUS only)."),
        "geoglows_vpu": ConfigOption(None, "text", "GEOGLOWS VPU (e.g. 714) for use_warning_flags_to_download_dem."),
    },
    "Streams": {
        "flowline": ConfigOption(None, "path", "Stream network file (shapefile, GeoPackage, GeoParquet, ...)."),
        "source_flowlines": ConfigOption([], "path_list", "Stream network files to take the domain's streams from. "
                                                          "Only those overlapping the DEM are read."),
        "stream_id_field": ConfigOption(None, "text", "Stream attribute with each reach's ID. Default: LINKNO for "
                                                      "GEOGLOWS, COMID for NWM; required for CUSTOM."),
        "downstream_id_field": ConfigOption(None, "text", "Stream attribute with the ID of the reach downstream. "
                                                          "Default: DSLINKNO for GEOGLOWS, TOCOMID for NWM; required "
                                                          "for CUSTOM."),
        "StrmOrder_Field": ConfigOption(None, "text", "Stream attribute with each reach's stream order."),
        "StrmOrder_Lower": ConfigOption(None, "int", "Map only reaches of at least this stream order. Not used by "
                                                     "Curve2Flood-FLDPLNpy."),
        "StrmOrder_Upper": ConfigOption(None, "int", "Map only reaches of at most this stream order. Not used by "
                                                     "Curve2Flood-FLDPLNpy."),
        "exclude": ConfigOption([], "int_list", "IDs of reaches to leave out."),
        "drop_multilinestrings": ConfigOption(False, "bool", "Leave out reaches whose line is in several pieces."),
        "lake_filter_json": ConfigOption(None, "path", "GEOGLOWS JSON of the reaches inside lakes. Those reaches are "
                                                       "left out."),
        "lakes": ConfigOption(None, "path", "Lake polygons. Lakes are kept out of burned channels and cut the moved "
                                            "stream network. Put the sea in ocean."),
        "ocean": ConfigOption(None, "path", "Ocean polygons, which the burned DEM has no data in (burn_streams): the "
                                            "streams end where they first reach the sea, and the reaches below are "
                                            "left out. A DEM's cells at 0 m aren't taken for the sea."),
        "move_stream_network_to_thalweg": ConfigOption(False, "bool", "Move each reach onto the DEM's flow path. "
                                                                      "Required by Curve2Flood-FLDPLNpy."),
        "new_strm_threshold_km2": ConfigOption(1, "float", "Smallest drainage area (km2) of the DEM flow paths that "
                                                           "reaches are moved onto."),
        "burn_streams": ConfigOption(False, "bool", "Burn the streams into the DEM before mapping."),
        "relocate_misplaced_streams": ConfigOption(True, "bool", "When moving the streams, first move reaches "
                                                                "the DEM shows are off their river onto the flow "
                                                                "paths their lines drain along."),
        "use_dem_derived_channel_mask": ConfigOption(False, "bool", "Find the channel the DEM shows along the "
                                                                    "streams: when burning, smooth the water it "
                                                                    "flattens and leave the lines unlowered on it, "
                                                                    "and add the channel to the water mask."),
        "streams_as_parquet": ConfigOption(False, "bool", "Write stream networks as GeoParquet instead of GeoPackage."),
        "minimize_output_files": ConfigOption(False, "bool", "Delete the domain's source and relocated streams "
                                                             "once the moved stream network is written."),
    },
    "Flows": {
        "streamflow_source": ConfigOption(str(StreamflowSource.GEOGLOWS_V2), "text", "Where the flows come from: "
                                                                                     "GEOGLOWS V2 (GEOGLOWS also "
                                                                                     "means V2), an NWM forecast, or "
                                                                                     "CUSTOM: only your files, with "
                                                                                     "the ID fields given.",
                                          tuple(str(source) for source in StreamflowSource)),
        "nwm_api_key": ConfigOption(None, "text", "NWM API key. Required for the NWM sources."),
        "floodmap_mode": ConfigOption(FloodMapMode.FORECAST.value, "text", "Which flows to map: a forecast, "
                                                                           "user_flow_files, or return_periods.",
                                      tuple(mode.value for mode in FloodMapMode)),
        "user_flow_files": ConfigOption([], "path_list", "CSV flow files to map in user mode: reach ID in the first "
                                                         "column, flows (m3/s) in the others."),
        "return_periods": ConfigOption([], "int_list", "Return periods to map in return_period mode: any of 2, 5, 10, "
                                                       "25, 50 and 100. GEOGLOWS or CUSTOM."),
        "reanalysis_file": ConfigOption(None, "path", "Flow statistics file (a COMID column, then flow columns) to use "
                                                      "instead of building one."),
        "return_period_file": ConfigOption(None, "text", "Return period flows to use instead of GEOGLOWS' (a CUSTOM "
                                                         "source's, unless reanalysis_file is given): a local path or "
                                                         "URL to a CSV, Parquet, Zarr or NetCDF file."),
        "return_period_variables": ConfigOption(None, "str_list", "Variables of the return period dataset to read; "
                                                                  "the highest is used. Default: all of them."),
        "include_fdc": ConfigOption(True, "bool", "Add flow duration curve flows (p_exceed_*) to the flow statistics. "
                                                  "GEOGLOWS, or CUSTOM with fdc_file."),
        "fdc_file": ConfigOption(None, "text", "Flow duration curves to use instead of GEOGLOWS', read like "
                                               "return_period_file."),
        "fdc_variable": ConfigOption(None, "text", "Variable of the flow duration curve dataset to read. Default: "
                                                   "hourly_annual, or its only one."),
        "reanalysis_storage_options": ConfigOption(None, "dict", "fsspec options for reading a remote "
                                                                 "return_period_file or fdc_file. Default: "
                                                                 "anonymous S3."),
        "raise_errors_if_river_ids_missing": ConfigOption(True, "bool", "Stop if reaches are missing from the flow "
                                                                        "datasets. If false, leave them out."),
        "specified_bathyflow_field": ConfigOption("rp2", "text", "Flow statistic that sets the baseflow the "
                                                                 "bathymetry is estimated for.",
                                                  tuple(FLOW_FIELD_OPTIONS)),
        "specified_highflow_field": ConfigOption("rp100_premium", "text", "Flow statistic of the highest flow the "
                                                                          "rating curves cover.",
                                                 tuple(FLOW_FIELD_OPTIONS)),
        "q_baseflow_threshold": ConfigOption(None, "float", "Leave out reaches whose baseflow is below this (m3/s)."),
    },
    "Forecasts": {
        "forensic_forecast_date": ConfigOption(None, "text", "Map the forecast issued on this date (YYYYMMDD, UTC). "
                                                             "Default: the newest forecast."),
        "forensic_forecast_hour": ConfigOption(None, "text", "Hour (UTC) of a past NWM forecast: 00-23 for short "
                                                             "range, 00, 06, 12 or 18 for medium range, 00 for long "
                                                             "range."),
        "remove_old_forecast_files": ConfigOption(False, "bool", "Delete this DEM's forecast flood maps once they are "
                                                                 "age_of_forecast_days old."),
        "age_of_forecast_days": ConfigOption(7, "int", "Age in days at which remove_old_forecast_files deletes a "
                                                       "forecast flood map."),
    },
    "Land cover and roughness": {
        "land_cover_cache": ConfigOption([], "path_list", "Local ESA WorldCover tiles to use instead of downloading "
                                                          "them, named by tile, e.g. N36W120.tif."),
        "land_watervalue": ConfigOption(80, "int", "Land cover value of open water (80 in ESA WorldCover)."),
        "mannings_text_file": ConfigOption(None, "path", "Manning's n of each land cover class. Default: a table for "
                                                         "ESA WorldCover."),
    },
    "Bathymetry": {
        "disable_bathymetry": ConfigOption(False, "bool", "Skip bathymetry and map floods on the DEM as it is."),
        "clean_dem": ConfigOption(False, "bool", "Smooth the stream channels with the DEM cleaner first."),
        "use_specified_depth_for_bathy_mask": ConfigOption(True, "bool", "Outline the DEM cleaner's channels by "
                                                                         "flooding the streams to set depths "
                                                                         "(specify_depths_for_bathy_mask) instead of "
                                                                         "using the water mask, and have ARC find the "
                                                                         "banks without the water mask."),
        "specify_depths_for_bathy_mask": ConfigOption([0.1], "float_list", "Flood depths (m): one, or two with "
                                                                           "clean_dem (the DEM cleaner's outline, then "
                                                                           "the bathymetry's). Only the DEM cleaner's "
                                                                           "has an effect: ARC's channel is burned in "
                                                                           "clipped to no outline."),
        "bathy_use_banks": ConfigOption(False, "bool", "Estimate bathymetry from the banks' elevation instead of the "
                                                       "water surface."),
        "find_banks_based_on_landcover": ConfigOption(True, "bool", "Put the banks where water land cover ends, for "
                                                                    "streams in water."),
        "use_power_laws_for_bathymetry": ConfigOption(False, "bool", "Size the channel from its drainage area with "
                                                                     "power laws instead of from the baseflow."),
        "coefficient_depth": ConfigOption(0.27, "float", "Power law: depth (m) = coefficient_depth * area "
                                                         "(km2)^exponent_depth."),
        "exponent_depth": ConfigOption(0.21, "float", "Power law exponent for depth."),
        "coefficient_width": ConfigOption(2.44, "float", "Power law: width (m) = coefficient_width * area "
                                                         "(km2)^exponent_width."),
        "exponent_width": ConfigOption(0.34, "float", "Power law exponent for width."),
        "area_km2_field": ConfigOption("DSContArea_km2", "text", "Stream attribute with each reach's drainage area in "
                                                                 "km2."),
        "area_m2_field": ConfigOption("DSContArea", "text", "Stream attribute with drainage area in m2, used when "
                                                            "area_km2_field is missing."),
        "slope_low_percentile": ConfigOption(25, "int", "Raise a stream cell's slope to at least this percentile of "
                                                        "its reach's slopes."),
        "slope_high_percentile": ConfigOption(75, "int", "Lower a stream cell's slope to at most this percentile of "
                                                         "its reach's slopes."),
        "create_reach_average_curve_file": ConfigOption(False, "bool", "Write one rating curve per reach instead of "
                                                                       "one per stream cell."),
        "standard_step": ConfigOption(False, "bool", "Link ARC's water surfaces along the stream network by the "
                                                     "standard-step method, instead of each cross section's normal "
                                                     "depth."),
        "valley_sections": ConfigOption(False, "bool", "Take ARC's overbank flow on lines out from each cross "
                                                       "section's banks square to the valley, at the valley's slope: "
                                                       "over a meandering river's floodplain, the water runs down the "
                                                       "valley, not along the channel."),
        "valley_dir_widths": ConfigOption(16.0, "float", "How far along the stream, in channel widths either way, ARC "
                                                         "takes the valley's direction and the stream's sinuosity "
                                                         "over, for valley_sections: about a meander's length, which "
                                                         "is 10 to 14 widths. 0 uses valley_dir_dist everywhere."),
        "valley_dir_dist": ConfigOption(40, "int", "How far along the stream, in cells either way, ARC takes the "
                                                   "valley's direction and the stream's sinuosity over, for "
                                                   "valley_sections, where a reach has no banks to give its width "
                                                   "(or with valley_dir_widths 0). At least 2."),
        "bathy_args": ConfigOption(DEFAULT_BATHY_ARGS, "dict", "More ARC settings for the bathymetry run. Replaces "
                                                               "the whole default set."),
    },
    "Outputs": {
        "make_vdt": ConfigOption(True, "bool", "Write ARC's VDT database (velocity, depth and top width per stream "
                                               "cell and flow). The flood maps use it, or else the curve file; "
                                               "Curve2Flood-FLDPLNpy needs it."),
        "make_curvefile": ConfigOption(True, "bool", "Write ARC's curve file (rating curve per stream cell). FIST "
                                                     "inputs are made from it."),
        "make_ap_database": ConfigOption(True, "bool", "Write ARC's area-perimeter database."),
        "make_cross_section_file": ConfigOption(False, "bool", "Write ARC's cross sections."),
        "make_representative_cross_section_file": ConfigOption(False, "bool", "Write a representative cross section "
                                                                              "per reach, for routing models."),
        "vdt_file_extension": ConfigOption("txt", "text", "File type of ARC's VDT and AP databases.",
                                           ("txt", "csv", "parquet")),
        "make_depth_maps": ConfigOption(True, "bool", "Write flood depth maps."),
        "make_velocity_maps": ConfigOption(True, "bool", "Write flood velocity maps."),
        "make_wse_maps": ConfigOption(True, "bool", "Write water surface elevation maps."),
        "flood_waterlc_and_strm_cells": ConfigOption(False, "bool", "Show water land cover and stream cells as "
                                                                    "flooded."),
        "flood_bathy_mask": ConfigOption(False, "bool", "Show the bathymetry water mask's cells as flooded."),
        "floodmap_identifier": ConfigOption("", "text", "Text added to the flood map file names."),
        "floodmap_args": ConfigOption(DEFAULT_FLOODMAP_ARGS, "dict", "More Curve2Flood settings for the flood maps. "
                                                                     "Replaces the whole default set."),
        "overwrite_floodmaps": ConfigOption(True, "bool", "Remake flood maps that already exist."),
        "make_fist_inputs": ConfigOption(True, "bool", "Write FIST inputs (GeoJSON water surface points) for each "
                                                       "forecast or user flow file. Needs make_curvefile."),
        "estimate_consequences": ConfigOption(False, "bool", "Estimate flood damage with Go-Consequences (needs "
                                                             "Docker; CONUS only). Turns on make_depth_maps."),
        "folder_paths": ConfigOption(None, "dict", "Where each output folder goes, as path templates. Default: "
                                                   "output_dir/name/<folder>."),
        "short_file_names": ConfigOption(False, "bool", "Drop the <source>_<DEM>_ prefix from file names."),
        "file_names": ConfigOption(None, "dict", "New names for kinds of output files."),
        "use_parquet": ConfigOption(False, "bool", "Write the FLDPLN stream info table as Parquet instead of CSV."),
        "use_yaml": ConfigOption(False, "bool", "Write ARC and Curve2Flood input files as YAML instead of text."),
        "compression": ConfigOption("LZW", "text", "Compression of the rasters nencarta writes, e.g. LZW, DEFLATE or "
                                                   "ZSTD."),
    },
    "Flood mapper": {
        "mapper": ConfigOption(Mapper.CURVE2FLOOD_KERNEL_WEIGHTED.value, "text", "Flood mapper. FloodSpreader needs "
                                                                                 "floodspreader.py in the nencarta "
                                                                                 "package folder.",
                               tuple(Mapper.list_names())),
        "fldpln_dh": ConfigOption(0.5, "float", "FLDPLN depth step (m). Used by the iterative solver."),
        "fldpln_min_depth": ConfigOption(0.1, "float", "FLDPLN: least depth (m) the iterative solver counts a flooded "
                                                       "cell as when it decides where a spill stops (FLDPLN's fldmn). "
                                                       "The library stores each depth to flood as found, however "
                                                       "small."),
        "fldpln_min_flood_depth": ConfigOption(0.01, "float", "FLDPLN: shallowest water (m) that floods a cell on the "
                                                              "flood maps (Curve2Flood's FLDPLN_Min_Flood_Depth). "
                                                              "0 floods every cell the water reaches."),
        "fldpln_max_depth": ConfigOption(25.0, "float", "FLDPLN: deepest flood (m) the library covers."),
        "fldpln_spill_decay": ConfigOption(0.0, "float", "FLDPLN: head (m) water spilling over a divide loses per "
                                                         "cell down its path, so a spill on flat ground ends where "
                                                         "its water runs out. 0 for none; 0.005 is a good starting "
                                                         "point. The exact solver takes it off on every step down, "
                                                         "whichever neighbour that is, as a spill loses it down its "
                                                         "path."),
        "fldpln_max_cells_per_segment": ConfigOption(None, "int", "FLDPLN: most floodplain cells one stream segment "
                                                                  "may flood. Saves memory on flat ground. Exact "
                                                                  "solver only."),
        "fldpln_keep_spilling": ConfigOption(False, "bool", "FLDPLN: keep routing water that spills over low divides. "
                                                            "Suits DEMs of 10 m or finer."),
        "fldpln_solver": ConfigOption("iterative", "text", "FLDPLN solver: iterative (steps up by fldpln_dh and "
                                                           "spills once per step, as FLDPLN advises for 30 m DEMs) "
                                                           "or exact (one sweep; spills as far as the water can run "
                                                           "downhill, and ignores fldpln_dh).",
                                      ("iterative", "exact")),
        "fldpln_parallel": ConfigOption(False, "bool", "FLDPLN: build the library with several processes. Cannot be "
                                                       "used with parallel."),
    },
    "Running": {
        "overwrite": ConfigOption(False, "bool", "Remake intermediate files that already exist."),
        "parallel": ConfigOption(False, "bool", "Process the DEMs in parallel."),
        "num_workers": ConfigOption(None, "int", "Number of processes when parallel."),
        "profile": ConfigOption(True, "bool", "Print how long each step took (serial runs only)."),
        "quiet": ConfigOption(False, "bool", "Hide ARC and Curve2Flood progress output."),
        "raise_errors_if_nothing_in_domain": ConfigOption(True, "bool", "Stop with an error when a DEM has no streams "
                                                                        "or data. If false, skip it."),
    },
}

CONFIG_OPTIONS = {key: option for options in CONFIG_SECTIONS.values() for key, option in options.items()}

DEFAULT_CONFIG = {key: option.default for key, option in CONFIG_OPTIONS.items()}
