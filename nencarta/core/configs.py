import os
from copy import deepcopy
from datetime import datetime
from collections.abc import Sequence

from nencarta.logger import LOG
from nencarta.core.defaults import DEFAULT_CONFIG, VALID_RETURN_PERIODS
from nencarta.core.enumerations import FloodMapMode, Mapper, StreamflowSource

# The reach ID and downstream reach ID attributes of each source's stream networks; a custom network's must be given
STREAM_ID_FIELDS = {"GEOGLOWS": ("LINKNO", "DSLINKNO"), "NWM": ("COMID", "TOCOMID"), "CUSTOM": (None, None)}

# The hours each NWM forecast is issued at
NWM_FORECAST_HOURS = {
    StreamflowSource.NWM_SHORT_RANGE: tuple(f"{hour:02d}" for hour in range(24)),
    StreamflowSource.NWM_MEDIUM_RANGE: ("00", "06", "12", "18"),
    StreamflowSource.NWM_LONG_RANGE: ("00",),
}

def strip_quotes(text):
    """``text`` without the quotes around it, as a path copied from a file browser has; other values as they are."""
    if not isinstance(text, str):
        return text
    text = text.strip()
    while len(text) >= 2 and text[0] == text[-1] and text[0] in "'\"":
        text = text[1:-1].strip()
    return text

def norm_or_none(path: str):
    path = strip_quotes(path)
    return os.path.normpath(path) if path else None


class NencartaConfig:
    """
    One watershed's settings: a config dictionary (keys in core/defaults.py), checked, with defaults filled in for
    missing keys and keys set to None.
    """

    def __init__(self, config_dict: dict):
        self.config_dict = config_dict
        self.setup()

    def validate_forecast_hour(self):
        """Return forensic_forecast_hour as two digits, e.g. '06', or None."""
        hour = self.get("forensic_forecast_hour")
        if hour is None or hour == "":
            return None
        try:
            hour = int(hour)
        except (TypeError, ValueError):
            raise ValueError(f"Invalid forensic_forecast_hour: {hour}")
        if hour not in range(24):
            raise ValueError(f"Invalid forensic_forecast_hour: {hour}. It must be between 0 and 23.")
        return f"{hour:02d}"

    def validate_forecast_hours(self):
        """Check forensic_forecast_hour is an hour the selected NWM forecast is issued at."""
        allowed = NWM_FORECAST_HOURS.get(self.streamflow_source)
        if self.forensic_forecast_hour and allowed and self.forensic_forecast_hour not in allowed:
            raise ValueError(
                f"Watershed '{self.watershed_name}' requires 'forensic_forecast_hour' to be one of "
                f"{', '.join(allowed)} when 'streamflow_source' is '{self.streamflow_source}'."
            )

    def validate_nwm_api_key(self):
        """Require an API key for NWM sources."""
        if self.streamflow_source.is_nwm() and not self.nwm_api_key:
            raise ValueError(f"Watershed '{self.watershed_name}' requires 'nwm_api_key' (or the NWM_API_KEY "
                             "environment variable) when 'streamflow_source' is NWM.")

    def validate_custom_source(self):
        """
        A custom source downloads nothing: it needs the stream network's ID fields, and the flow statistics from
        reanalysis_file, or return_period_file (with fdc_file for p_exceed flows). Its maps are made from user flow
        files or return periods; there are no forecasts.
        """
        if self.streamflow_source != StreamflowSource.CUSTOM:
            return
        name = f"Watershed '{self.watershed_name}'"
        if not self.stream_id_field or not self.downstream_id_field:
            raise ValueError(f"{name}: a custom 'streamflow_source' needs 'stream_id_field' and 'downstream_id_field'.")
        if self.floodmap_mode == FloodMapMode.FORECAST:
            raise ValueError(f"{name}: forecasts come from GEOGLOWS or NWM; with a custom 'streamflow_source', map "
                             "'user_flow_files' or 'return_periods' instead.")
        if self.reanalysis_file:
            return
        if not self.return_period_file:
            raise ValueError(f"{name}: a custom 'streamflow_source' takes its flow statistics from 'reanalysis_file', "
                             "or from 'return_period_file' (and 'fdc_file' for p_exceed flows).")
        fdc_fields = [field for field in (self.specified_bathyflow_field, self.specified_highflow_field)
                      if str(field).startswith("p_exceed")]
        if fdc_fields and not (self.fdc_file and self.include_fdc):
            raise ValueError(f"{name}: {', '.join(fdc_fields)} come from a flow duration curve, so a custom "
                             "'streamflow_source' needs 'fdc_file' (and 'include_fdc') for them.")

    def validate_specified_depths(self):
        """Check there are as many bathymetry mask depths as the run uses: two with clean_dem, otherwise one."""
        if self.use_specified_depth_for_bathy_mask:
            if not self.specify_depths_for_bathy_mask or not isinstance(self.specify_depths_for_bathy_mask, list) or len(self.specify_depths_for_bathy_mask) not in {1, 2}:
                raise ValueError(f"Watershed '{self.watershed_name}' requires 'specify_depths_for_bathy_mask' as a list of 1-2 floats when 'use_specified_depth_for_bathy_mask' is True.")
            elif len(self.specify_depths_for_bathy_mask) != 2 and self.clean_dem:
                raise ValueError(f"Watershed '{self.watershed_name}' requires 'specify_depths_for_bathy_mask' as a list of two floats when 'clean_dem' is True.")
            elif len(self.specify_depths_for_bathy_mask) != 1 and not self.clean_dem:
                raise ValueError(f"Watershed '{self.watershed_name}' requires 'specify_depths_for_bathy_mask' as a list of one float when 'clean_dem' is False.")

    def validate_forecast_date(self):
        """Return forensic_forecast_date (YYYYMMDD) or None."""
        date = self.get("forensic_forecast_date")
        if not date:
            return None
        for date_format in ('%Y%m%d', '%Y-%m-%d %H:%M:%S %Z'):
            try:
                datetime.strptime(date, date_format)
                return date
            except ValueError:
                pass
        raise ValueError(f"Invalid forensic_forecast_date format: {date}")

    def verify_required_keys(self):
        """Check for a name, an output folder, a DEM source and a stream network."""
        required_keys = ["name", "output_dir",]
        for key in required_keys:
            if key not in self.config_dict:
                raise KeyError(f"Missing required key in watershed {self.get('name', 'unknown')}: {key}")
        if not any(key in self.config_dict for key in ["dem", "source_dems", "dem_dir"]):
            raise KeyError(f"Watershed '{self.get('name', 'unknown')}' must have at least one of the following keys: 'dem', 'source_dems', or 'dem_dir'.")
        if not any(key in self.config_dict for key in ["flowline", "source_flowlines"]):
            raise KeyError(f"Watershed '{self.get('name', 'unknown')}' must have at least one of the following keys: 'flowline', or 'source_flowlines'.")

    def warn_about_unknown_keys(self):
        """Warn about keys nencarta doesn't use, such as misspelt or retired ones."""
        unknown = sorted(set(self.config_dict) - set(DEFAULT_CONFIG))
        if unknown:
            LOG.warning(f"Watershed '{self.watershed_name}': ignoring unknown keys {unknown}.")

    def warn_if_multiple_primary_dems(self):
        """Warn when both dem and dem_dir are given."""
        if self.get_path("dem") and self.get_path("dem_dir"):
            LOG.warning(
                f"Watershed '{self.watershed_name}' has both 'dem' and 'dem_dir'. "
                "'dem_dir' takes precedence in process_watershed; provide only one "
                "to avoid ambiguity."
            )

    def validate_user_floodmaps(self):
        """Return the flood map mode and the user flow files, as a list of paths."""
        floodmap_mode = self.get("floodmap_mode", "forecast")
        floodmap_mode = FloodMapMode.from_string(floodmap_mode)

        user_flow_files = self.get("user_flow_files", [])
        if isinstance(user_flow_files, str):
            user_flow_files = [user_flow_files]

        if floodmap_mode == FloodMapMode.USER and not isinstance(user_flow_files, list):
            raise ValueError(f"Watershed '{self.get('name', 'unknown')}' requires 'user_flow_files' as either a filepath string or a list of file paths when 'floodmap_mode' is 'user'.")

        if user_flow_files:
            user_flow_files = [os.path.normpath(f) for f in user_flow_files]

        return floodmap_mode, user_flow_files

    def validate_bbox(self):
        """Return bbox as four floats, or None."""
        bbox = self.get("bbox")
        if bbox:
            if not isinstance(bbox, Sequence) or len(bbox) != 4:
                raise ValueError(f"Watershed '{self.get('name', 'unknown')}' has invalid 'bbox': {bbox}. Must be a sequence of four floats [minx, miny, maxx, maxy].")
            try:
                bbox = [float(coord) for coord in bbox]
            except ValueError:
                raise ValueError(f"Watershed '{self.get('name', 'unknown')}' has non-numeric values in 'bbox': {bbox}. All coordinates must be convertible to float.")
        return bbox

    def validate_return_periods(self):
        """Return the return periods as ints, checked to be ones nencarta has flows for."""
        return_periods = self.get("return_periods")
        if return_periods and self.floodmap_mode == FloodMapMode.RETURN_PERIOD:
            if isinstance(return_periods, str) or not isinstance(return_periods, Sequence):
                return_periods = [return_periods]
            try:
                return_periods = [int(rp) for rp in return_periods]
            except ValueError:
                raise ValueError(f"Watershed '{self.get('name', 'unknown')}' has non-numeric values in 'return_periods': {return_periods}. All return periods must be whole numbers.")

            invalid = [rp for rp in return_periods if rp not in VALID_RETURN_PERIODS]
            if invalid:
                raise ValueError(f"Invalid return periods {invalid}. Must be a subset of {set(VALID_RETURN_PERIODS)}.")

        return return_periods

    def get(self, key: str, default=None):
        """Return a copy of a config value, or its default when it is missing or None."""
        if default is None and key in DEFAULT_CONFIG:
            default = DEFAULT_CONFIG[key]
        value = self.config_dict.get(key, default)
        if value is None and default is not None:
            return deepcopy(default)
        return deepcopy(value)

    def get_path(self, key: str, default=None):
        """Return a normalized path string or None."""
        value = self.get(key, default)
        return norm_or_none(value)

    def get_num_workers(self):
        """Return num_workers as a positive int, or None."""
        num_workers = self.get("num_workers")
        if num_workers is None:
            return None
        try:
            num_workers = int(num_workers)
            if num_workers < 1:
                raise ValueError
        except (TypeError, ValueError):
            raise ValueError(f"Invalid 'num_workers' value: {num_workers}. Must be a positive integer.")
        return num_workers

    def setup(self):
        """Check the config and set its values as attributes."""
        self.verify_required_keys()
        self.watershed_name: str = self.get("name")
        self.warn_about_unknown_keys()
        self.warn_if_multiple_primary_dems()

        self.streamflow_source: StreamflowSource = StreamflowSource.from_string(self.get("streamflow_source"))
        source_id_fields = STREAM_ID_FIELDS["NWM" if self.streamflow_source.is_nwm() else
                                            "CUSTOM" if self.streamflow_source == StreamflowSource.CUSTOM else "GEOGLOWS"]
        self.stream_id_field: str = self.get("stream_id_field") or source_id_fields[0]
        self.downstream_id_field: str = self.get("downstream_id_field") or source_id_fields[1]
        self.forensic_forecast_date: str | None = self.validate_forecast_date()
        self.forensic_forecast_hour: str | None = self.validate_forecast_hour()
        self.validate_forecast_hours()

        self.nwm_api_key: str = self.get("nwm_api_key") or os.getenv("NWM_API_KEY")
        self.validate_nwm_api_key()

        self.use_specified_depth_for_bathy_mask: bool = self.get("use_specified_depth_for_bathy_mask")
        self.specify_depths_for_bathy_mask: list[float] = self.get("specify_depths_for_bathy_mask")
        self.clean_dem: bool = self.get("clean_dem")
        self.validate_specified_depths()

        self.floodmap_mode, self.user_flow_files = self.validate_user_floodmaps()
        self.return_periods = self.validate_return_periods()

        if not self.get("make_depth_maps") and self.get('estimate_consequences'):
            LOG.warning(f"Watershed '{self.get('name', 'unknown')}': 'make_depth_maps' is False but 'estimate_consequences' is True. Setting 'make_depth_maps' to True.")
            self.config_dict['make_depth_maps'] = True

        self.dem_filter: str = self.get("dem_filter")
        if not self.dem_filter:
            self.dem_filter = "*"

        self.move_stream_network_to_thalweg: bool = self.get("move_stream_network_to_thalweg")
        try:
            self.new_strm_threshold_km2: float = float(self.get("new_strm_threshold_km2"))
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Watershed '{self.watershed_name}': 'new_strm_threshold_km2' must be a number, "
                             f"not {self.get('new_strm_threshold_km2')!r}.") from exc

        self.mapper: Mapper = Mapper.from_string(self.get('mapper'))

        self.streams_as_parquet: bool = self.get("streams_as_parquet")

        # Flood map names end with '_<floodmap_identifier>'
        self.floodmap_identifier: str = self.get("floodmap_identifier")
        self.floodmap_id: str = self.floodmap_identifier
        if self.floodmap_id:
            self.floodmap_id = f"_{self.floodmap_id}" if not self.floodmap_id.startswith("_") else self.floodmap_id

        self.bbox: tuple[float, ...] = self.validate_bbox()
        self.flowline: str = self.get_path("flowline")
        self.source_flowlines: list[str] = self.get("source_flowlines")
        self.dem: str = self.get_path("dem")
        self.dem_dir: str = self.get_path("dem_dir")
        self.output_dir: str = self.get_path("output_dir")
        self.folder_paths: dict[str, str] | None = self.get("folder_paths")
        if self.folder_paths is not None and not isinstance(self.folder_paths, dict):
            raise ValueError(f"Watershed '{self.get('name', 'unknown')}': 'folder_paths' must be a dictionary of folder name to path template.")
        self.short_file_names: bool = self.get("short_file_names")
        self.file_names: dict[str, str] | None = self.get("file_names")
        if self.file_names is not None and not isinstance(self.file_names, dict):
            raise ValueError(f"Watershed '{self.get('name', 'unknown')}': 'file_names' must be a dictionary of file to name.")
        self.disable_bathymetry: bool = self.get("disable_bathymetry")
        self.bathy_use_banks: bool = self.get("bathy_use_banks")
        self.flood_waterlc_and_strm_cells: bool = self.get("flood_waterlc_and_strm_cells")
        self.flood_bathy_mask: bool = self.get("flood_bathy_mask")
        self.land_watervalue: int = self.get("land_watervalue")
        self.overwrite: bool = self.get("overwrite")
        self.age_of_forecast_days: int = self.get("age_of_forecast_days")
        self.find_banks_based_on_landcover: bool = self.get("find_banks_based_on_landcover")
        self.create_reach_average_curve_file: bool = self.get("create_reach_average_curve_file")
        self.standard_step: bool = self.get("standard_step")
        self.use_warning_flags_to_download_dem: bool = self.get("use_warning_flags_to_download_dem")
        self.geoglows_vpu: str = self.get("geoglows_vpu")
        self.specified_bathyflow_field: str = self.get("specified_bathyflow_field")
        self.specified_highflow_field: str = self.get("specified_highflow_field")
        self.StrmOrder_Field: str = self.get("StrmOrder_Field")
        self.StrmOrder_Lower: int = self.get("StrmOrder_Lower")
        self.StrmOrder_Upper: int = self.get("StrmOrder_Upper")
        self.q_baseflow_threshold: float = self.get("q_baseflow_threshold")
        self.lake_filter_json: str = self.get_path("lake_filter_json")
        self.lakes: str = self.get_path("lakes")
        self.estimate_consequences: bool = self.get("estimate_consequences")
        self.overwrite_floodmaps: bool = self.get("overwrite_floodmaps")
        self.remove_old_forecast_files: bool = self.get("remove_old_forecast_files")
        self.make_fist_inputs: bool = self.get("make_fist_inputs")
        self.make_vdt: bool = self.get("make_vdt")
        self.make_curvefile: bool = self.get("make_curvefile")
        self.make_ap_database: bool = self.get("make_ap_database")
        self.make_cross_section_file: bool = self.get("make_cross_section_file")
        self.make_representative_cross_section_file: bool = self.get("make_representative_cross_section_file")
        self.vdt_file_extension: str = self.get("vdt_file_extension")
        self.mannings_text_file: str = self.get("mannings_text_file")
        self.bathy_args: dict = self.get("bathy_args")
        self.floodmap_args: dict = self.get("floodmap_args")
        self.make_depth_maps: bool = self.get("make_depth_maps")
        self.make_velocity_maps: bool = self.get("make_velocity_maps")
        self.make_wse_maps: bool = self.get("make_wse_maps")
        self.quiet: bool = self.get("quiet")
        self.source_dems: list = self.get("source_dems")
        self.buffer: bool = self.get("buffer")
        self.buffer_distance: float = self.get("buffer_distance")
        self.use_vrt: bool = self.get("use_vrt")
        self.raise_errors_if_nothing_in_domain: bool = self.get("raise_errors_if_nothing_in_domain")
        self.land_cover_cache: list[str] = self.get("land_cover_cache")
        self.use_parquet: bool = self.get("use_parquet")
        self.use_yaml: bool = self.get("use_yaml")
        # The date and hour of the newest forecast, once it has been found
        self.forecastdate: str = None
        self.forecasthour: str = None
        self.profile: bool = self.get("profile")
        self.parallel: bool = self.get("parallel")
        self.compression: str = self.get("compression")
        self.exclude: list[int] = self.get("exclude")
        self.reanalysis_file: str = self.get_path("reanalysis_file")
        self.return_period_file: str | None = strip_quotes(self.get("return_period_file"))
        self.return_period_variables: list[str] | None = self.get("return_period_variables")
        if isinstance(self.return_period_variables, str):
            self.return_period_variables = [self.return_period_variables]
        self.include_fdc: bool = self.get("include_fdc")
        self.fdc_file: str | None = strip_quotes(self.get("fdc_file"))
        self.fdc_variable: str | None = self.get("fdc_variable")
        self.reanalysis_storage_options: dict | None = self.get("reanalysis_storage_options")
        self.raise_errors_if_river_ids_missing: bool = self.get("raise_errors_if_river_ids_missing")
        self.fldpln_dh: float = self.get("fldpln_dh")
        self.fldpln_min_depth: float = self.get("fldpln_min_depth")
        self.fldpln_max_depth: float = self.get("fldpln_max_depth")
        self.fldpln_spill_decay: float = self.get("fldpln_spill_decay")
        self.fldpln_max_cells_per_segment: int | None = self.get("fldpln_max_cells_per_segment")
        self.fldpln_keep_spilling: bool = self.get("fldpln_keep_spilling")
        self.fldpln_parallel: bool = self.get("fldpln_parallel")
        self.fldpln_solver: str = self.get("fldpln_solver")
        self.project_to_utm: bool = self.get("project_to_utm")
        self.burn_streams: bool = self.get("burn_streams")
        self.relocate_misplaced_streams: bool = self.get("relocate_misplaced_streams")
        self.use_power_laws_for_bathymetry: bool = self.get("use_power_laws_for_bathymetry")
        self.area_m2_field: str = self.get("area_m2_field")
        self.area_km2_field: str = self.get("area_km2_field")
        self.num_workers = self.get_num_workers()
        self.slope_low_percentile: int = int(self.get("slope_low_percentile"))
        self.slope_high_percentile: int = int(self.get("slope_high_percentile"))
        self.coefficient_depth = self.get("coefficient_depth")
        self.coefficient_width = self.get("coefficient_width")
        self.exponent_depth = self.get("exponent_depth")
        self.exponent_width = self.get("exponent_width")
        self.use_dem_derived_channel_mask: bool = self.get("use_dem_derived_channel_mask")
        self.drop_multilinestrings: bool = self.get("drop_multilinestrings")
        self.minimize_output_files: bool = self.get("minimize_output_files")

        self.validate_custom_source()

        if self.parallel and self.fldpln_parallel:
            LOG.error("Both 'parallel' and 'fldpln_parallel' are set to True. Please set only one of these options to True.")
            raise ValueError("Both 'parallel' and 'fldpln_parallel' are set to True. Please set only one of these options to True.")

        if not self.make_vdt and not self.make_curvefile:
            raise ValueError(f"Watershed '{self.watershed_name}': the flood maps need ARC's VDT database or curve file, "
                             "so 'make_vdt' and 'make_curvefile' can't both be False.")
        if not self.make_vdt and self.mapper.is_curve2flood_fldpln_mapper():
            raise ValueError(f"Watershed '{self.watershed_name}': the Curve2Flood-FLDPLNpy mapper needs ARC's VDT "
                             "database, so 'make_vdt' must be True.")

        if self.use_power_laws_for_bathymetry and not all([self.coefficient_depth, self.coefficient_width, self.exponent_depth, self.exponent_width]):
            LOG.error("When 'use_power_laws_for_bathymetry' is True, all of 'coefficient_depth', 'coefficient_width', 'exponent_depth', and 'exponent_width' must be provided.")
            raise ValueError("When 'use_power_laws_for_bathymetry' is True, all of 'coefficient_depth', 'coefficient_width', 'exponent_depth', and 'exponent_width' must be provided.")

    def __repr__(self):
        return f"NencartaConfig(name={self.watershed_name}, output_dir={self.output_dir})"
