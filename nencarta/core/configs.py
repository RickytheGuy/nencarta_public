import os
from datetime import datetime
from collections.abc import Sequence

from nencarta.logger import LOG
from nencarta.core.enumerations import FloodMapMode, Mapper, StreamflowSource

def normalize_mapper_name(mapper: str | None) -> str:
    if mapper is None:
        return "Curve2Flood-Kernel Weighted"
    mapper = str(mapper).strip()
    if mapper == "":
        return "Curve2Flood-Kernel Weighted"
    return mapper


def is_curve2flood_fldpln_mapper(mapper: str | None) -> bool:
    return normalize_mapper_name(mapper) == "Curve2Flood-FLDPLNpy"

def norm_or_none(path: str):
    return os.path.normpath(path) if path else None

def float_or_none(value):
    if value in (None, ""):
        return None
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid q_baseflow_threshold: {value}") from exc

def get_streamids_from_source(streamflow_source: str):
    if streamflow_source.upper().startswith("NWM"):
        return 'COMID', 'TOCOMID'
    
    if streamflow_source.upper() == "GEOGLOWS":
        return 'LINKNO', 'DSLINKNO'
    
    LOG.error(f"streamflow_source {streamflow_source} not recognized, please use either 'NWM' or 'GEOGLOWS'")
    raise ValueError(f"streamflow_source {streamflow_source} not recognized, please use either 'NWM' or 'GEOGLOWS'")


class NencartaConfig:
    """Validated user configuration for a watershed run."""

    def __init__(self, config_dict: dict):
        self.config_dict = config_dict
        self.setup()

    
    def validate_forecast_hour(self):
        """Return a zero-padded forecast hour or None."""
        forensic_forecast_hour = self.get("forensic_forecast_hour")
        if forensic_forecast_hour:
            try:
                forensic_forecast_hour = int(forensic_forecast_hour)
                if forensic_forecast_hour not in range(0, 24):
                    raise ValueError("forensic_forecast_hour must be between 0 and 23")
                forensic_forecast_hour = f"{forensic_forecast_hour:02d}"
            except ValueError:
                raise ValueError(f"Invalid forensic_forecast_hour: {forensic_forecast_hour}")
                
        return forensic_forecast_hour

    def validate_forecast_hours(self):
        """Validate NWM forecast-hour rules for the selected source."""
        # if the streamflow_source is "NWM_short_range" the forensic_forecast_hour can be between 0 and 23, the forensic_forecast_hour must be provided as a two-digit string
        short_range_forecast_hours = [f"{i:02d}" for i in range(0, 24)]
        medium_range_forecast_hours = ["00", "06", "12", "18"]
        long_range_forecast_hours = ["00"]
        if self.streamflow_source == StreamflowSource.NWM_SHORT_RANGE and (self.forensic_forecast_hour and self.forensic_forecast_hour not in short_range_forecast_hours):
            raise ValueError(f"Watershed '{self.watershed_name}' requires 'forensic_forecast_hour' to be between 0 and 23 when 'streamflow_source' is 'NWM_short_range'.")
        # if the streamflow_source is "NWM_medium_range" the forensic_forecast_hour can be between one of 0, 6, 12, or 18
        if self.streamflow_source == StreamflowSource.NWM_MEDIUM_RANGE and (self.forensic_forecast_hour and self.forensic_forecast_hour not in medium_range_forecast_hours):
            raise ValueError(f"Watershed '{self.watershed_name}' requires 'forensic_forecast_hour' to be one of 0, 6, 12, or 18 when 'streamflow_source' is 'NWM_medium_range'.")
        # if the streamflow_source is "NWM_long_range" the forensic_forecast_hour can be only 0
        if self.streamflow_source == StreamflowSource.NWM_LONG_RANGE and (self.forensic_forecast_hour and self.forensic_forecast_hour not in long_range_forecast_hours):
            raise ValueError(f"Watershed '{self.watershed_name}' requires 'forensic_forecast_hour' to be 0 when 'streamflow_source' is 'NWM_long_range'.")

    def validate_nwm_api_key(self):
        """Require an API key for NWM sources."""
        if self.streamflow_source.upper().startswith("NWM") and not self.nwm_api_key:
            raise ValueError(f"Watershed '{self.watershed_name}' requires 'nwm_api_key' when 'streamflow_source' is NWM.")

    def validate_specified_depths(self):
        """Validate bathymetry mask depth count for cleaner and non-cleaner runs."""
        if self.use_specified_depth_for_bathy_mask:
            if not self.specify_depths_for_bathy_mask or not isinstance(self.specify_depths_for_bathy_mask, list) or len(self.specify_depths_for_bathy_mask) not in {1, 2}:
                raise ValueError(f"Watershed '{self.watershed_name}' requires 'specify_depths_for_bathy_mask' as a list of 1-2 floats when 'use_specified_depth_for_bathy_mask' is True.")
            elif len(self.specify_depths_for_bathy_mask) != 2 and self.clean_dem:
                raise ValueError(f"Watershed '{self.watershed_name}' requires 'specify_depths_for_bathy_mask' as a list of two floats when 'clean_dem' is True.")
            elif len(self.specify_depths_for_bathy_mask) != 1 and not self.clean_dem:
                raise ValueError(f"Watershed '{self.watershed_name}' requires 'specify_depths_for_bathy_mask' as a list of one float when 'clean_dem' is False.")

    def validate_forecast_date(self):
        """Return a valid forensic forecast date string or None."""
        # check if forensic_forecast_date and forensic_forecast_hour is provided in the watershed dictionary and if not set forensic_forecast_date=None
        # get forensic forecast date (string like "20231125" or "2023-11-25 06:00:00 UTC")
        forensic_forecast_date: str = self.get("forensic_forecast_date")
        if forensic_forecast_date:
            try:
                # First try YYYYMMDD format
                forensic_forecast_date_dt = datetime.strptime(forensic_forecast_date, '%Y%m%d')
            except ValueError:
                try:
                    # Fallback for full timestamp
                    forensic_forecast_date_dt = datetime.strptime(forensic_forecast_date, '%Y-%m-%d %H:%M:%S %Z')
                except ValueError:
                    raise ValueError(f"Invalid forensic_forecast_date format: {forensic_forecast_date}")
                
        else:
            forensic_forecast_date = None
            LOG.info("Forensic forecast date not provided; defaulting to None.")

        return forensic_forecast_date

    def verify_required_keys(self):
        """Validate the minimum inputs needed to build workspaces."""
        required_keys = ["name", "output_dir",]
        for key in required_keys:
            if key not in self.config_dict:
                raise KeyError(f"Missing required key in watershed {self.get('name', 'unknown')}: {key}")
        if not any(key in self.config_dict for key in ["dem", "source_dems", "dem_dir"]):
            raise KeyError(f"Watershed '{self.get('name', 'unknown')}' must have at least one of the following keys: 'dem', 'source_dems', or 'dem_dir'.")
        if not any(key in self.config_dict for key in ["flowline", "source_flowlines"]):
            raise KeyError(f"Watershed '{self.get('name', 'unknown')}' must have at least one of the following keys: 'flowline', or 'source_flowlines'.")

    def validate_user_floodmaps(self):
        """Validate and normalize user-supplied flood map flow files."""
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
        """Return a four-value bounding box as floats when provided."""
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
        """Validate requested return periods for return-period mode."""
        VALID_RETURN_PERIODS = {2, 5, 10, 25, 50, 100}
        return_periods = self.get("return_periods")
        if return_periods and self.floodmap_mode == FloodMapMode.RETURN_PERIOD:
            if not isinstance(return_periods, Sequence):
                raise ValueError(f"Watershed '{self.get('name', 'unknown')}' has invalid 'return_periods': {return_periods}. Must be a sequence of numeric values.")
            try:
                return_periods = [int(rp) for rp in return_periods]
            except ValueError:
                raise ValueError(f"Watershed '{self.get('name', 'unknown')}' has non-numeric values in 'return_periods': {return_periods}. All return periods must be convertible to float.")
        
            assert all(rp in VALID_RETURN_PERIODS for rp in return_periods), \
                f"Invalid return periods. Must be a subset of {VALID_RETURN_PERIODS}. Received: {return_periods}"
            
        return return_periods
        
    def get(self, key: str, default=None):
        """Return a config value, treating None as missing when a default exists."""
        value = self.config_dict.get(key, default)
        if value is None and default is not None:
            return default
        return value
    
    def get_path(self, key: str, default=None):
        """Return a normalized path string or None."""
        value = self.get(key, default)
        return norm_or_none(value)
        
    def setup(self):
        """Populate validated attributes used by the processing pipeline."""
        self.verify_required_keys()
        self.watershed_name: str = self.get("name")

        self.streamflow_source: StreamflowSource = StreamflowSource.from_string(self.get("streamflow_source", "GEOGLOWS"))
        self.stream_id_field, self.downstream_id_field = get_streamids_from_source(self.streamflow_source)
        self.forensic_forecast_date: str | None = self.validate_forecast_date()
        self.forensic_forecast_hour: int = self.validate_forecast_hour()
        self.validate_forecast_hours()

        self.nwm_api_key: str = self.get("nwm_api_key")
        self.validate_nwm_api_key()

        self.use_specified_depth_for_bathy_mask: bool = self.get("use_specified_depth_for_bathy_mask", True)
        self.specify_depths_for_bathy_mask: list[float] = self.get("specify_depths_for_bathy_mask")
        self.clean_dem: bool = self.get("clean_dem", False)
        self.validate_specified_depths()

        self.floodmap_mode, self.user_flow_files = self.validate_user_floodmaps()
        self.return_periods = self.validate_return_periods()

        if not self.get("make_depth_maps", True) and self.get('estimate_consequences', False):
            LOG.warning(f"Watershed '{self.get('name', 'unknown')}': 'make_depth_maps' is False but 'estimate_consequences' is True. Setting 'make_depth_maps' to True.")
            self.config_dict['make_depth_maps'] = True

        self.dem_filter: str = self.get("dem_filter", "*")
        if not self.dem_filter:
            self.dem_filter = "*"

        self.move_stream_network_to_thalweg: bool = self.get("move_stream_network_to_thalweg", False)
        self.stream_order_threshold: float = float_or_none(self.get("new_strm_threshold_km2", 5))  # TDX-Hydro uses 5 km2. It makes the conflation step easier if we can match
        if self.move_stream_network_to_thalweg and self.stream_order_threshold is None:
            raise ValueError(f"Watershed '{self.get('name', 'unknown')}': 'new_strm_threshold_km2' must be specified when moving stream network.")

        self.mapper: Mapper = Mapper.from_string(self.get('mapper'))
        
        if self.mapper.is_curve2flood_fldpln_mapper() and not self.move_stream_network_to_thalweg:
            LOG.warning(f"Move stream network was set to false, but the mapper is set to Curve2Flood-FLDPLNpy. Ensure that the input stream network is properly aligned with the DEM.")
        
        self.streams_as_parquet: bool = self.get("streams_as_parquet", False)
        if self.streams_as_parquet and self.mapper.is_curve2flood_fldpln_mapper():
            raise ValueError("The Curve2Flood-FLDPLNpy mapper is not currently compatible with stream networks stored as parquet files. Please set streams_as_parquet to false in order to proceed...")
        
        self.floodmap_id: str = self.get('floodmap_identifier', '')
        if self.floodmap_id:
            self.floodmap_id = f"_{self.floodmap_id}"if not self.floodmap_id.startswith("_") else self.floodmap_id
        else:
            self.floodmap_id = ''

        self.bbox: tuple[float, ...] = self.validate_bbox()
        self.flowline: str = self.get_path("flowline")
        self.source_flowlines: list[str] = self.get("source_flowlines", [])
        self.dem: str = self.get_path("dem")
        self.dem_dir: str = self.get_path("dem_dir")
        self.output_dir: str = self.get_path("output_dir")
        self.disable_bathymetry: bool = self.get("disable_bathymetry", False)
        self.bathy_use_banks: bool = self.get("bathy_use_banks", False)
        self.flood_waterlc_and_strm_cells: bool = self.get("flood_waterlc_and_strm_cells", False)
        self.land_watervalue: int = self.get("land_watervalue", 80)
        self.process_stream_network: bool = self.get("process_stream_network", False)
        self.age_of_forecast_days: int = self.get("age_of_forecast_days", 7)
        self.find_banks_based_on_landcover: bool = self.get("find_banks_based_on_landcover", True)
        self.create_reach_average_curve_file: bool = self.get("create_reach_average_curve_file", False)
        self.use_warning_flags_to_download_dem: bool = self.get("use_warning_flags_to_download_dem", False)
        self.geoglows_vpu: str = self.get("geoglows_vpu")
        self.specified_bathyflow_field: str = self.get("specified_bathyflow_field", "p_exceed_50")
        self.specified_highflow_field: str = self.get("specified_highflow_field", "rp100_premium")
        self.StrmOrder_Field: str = self.get("StrmOrder_Field")
        self.StrmOrder_Lower: int = self.get("StrmOrder_Lower")
        self.StrmOrder_Upper: int = self.get("StrmOrder_Upper")
        self.q_baseflow_threshold: float = self.get("q_baseflow_threshold")
        self.lake_filter_json: str = self.get_path("lake_filter_json")
        self.lakes: str = self.get_path("lakes")
        self.estimate_consequences: bool = self.get("estimate_consequences", False)
        self.overwrite_floodmaps: bool = self.get("overwrite_floodmaps", True)
        self.remove_old_forecast_files: bool = self.get("remove_old_forecast_files", False)
        self.make_fist_inputs: bool = self.get("make_fist_inputs", True)
        self.make_vdt: bool = self.get("make_vdt", True)
        self.make_curvefile: bool = self.get("make_curvefile", True)
        self.make_ap_database: bool = self.get("make_ap_database", True)
        self.make_cross_section_file: bool = self.get("make_cross_section_file", False)
        self.vdt_file_extension: str = self.get("vdt_file_extension", 'txt')
        self.mannings_text_file: str = self.get("mannings_text_file")
        self.bathy_args: dict = self.get("bathy_args", {})
        self.floodmap_args: dict = self.get("floodmap_args", {})
        self.make_depth_maps: bool = self.get("make_depth_maps", True)
        self.make_velocity_maps: bool = self.get("make_velocity_maps", True)
        self.make_wse_maps: bool = self.get("make_wse_maps", True)
        self.floodmap_identifier: str = self.get("floodmap_identifier", "")
        self.new_strm_threshold_km2: float = self.get("new_strm_threshold_km2", 5 ** 2) # TDX-Hydro uses 5 km2. It makes the conflation step easier if we can match
        self.min_match_score: float = self.get("min_match_score")
        self.quiet: bool = self.get("quiet", False)
        self.source_dems: list = self.get("source_dems", [])
        self.buffer: bool = self.get("buffer", False)
        self.buffer_distance: float = self.get("buffer_distance", 0.1)
        self.use_vrt: bool = self.get("use_vrt", False)
        self.raise_errors_if_nothing_in_domain: bool = self.get("raise_errors_if_nothing_in_domain", True)
        self.land_cover_cache: list[str] = self.get("land_cover_cache", [])
        self.use_parquet: bool = self.get("use_parquet", False)
        self.use_yaml: bool = self.get("use_yaml", False)
        self.forecastdate: str = None
        self.forecasthour: str = None
        self.profile: bool = self.get("profile", True)
        self.parallel: bool = self.get("parallel", False)
        self.compression: str = self.get("compression", "LZW")
        self.exclude: list[int] = self.get("exclude", [])
        self.reanalysis_file: str = self.get("reanalysis_file")
        self.fldpln_dh: float = self.get("fldpln_dh", 0.5)
        self.fldpln_min_depth: float = self.get("fldpln_min_depth", 0.1)
        self.fldpln_max_depth: float = self.get("fldpln_max_depth", 25.0)
        self.fldpln_keep_spilling: bool = self.get("fldpln_keep_spilling", False)
        self.fldpln_parallel: bool = self.get("fldpln_parallel", False)
        self.fldpln_max_wse_rise: float = self.get("fldpln_max_wse_rise", 0.01)
        self.project_to_utm: bool = self.get("project_to_utm", False)
        self.burn_streams: bool = self.get("burn_streams", False)
        self.use_power_laws_for_bathymetry: bool = self.get("use_power_laws_for_bathymetry", False)
        self.area_m2_field: str = self.get("area_m2_field", "DSContArea")
        self.area_km2_field: str = self.get("area_km2_field", "DSContArea_km2")

    def __repr__(self):
        return f"NencartaConfig(name={self.watershed_name}, output_dir={self.output_dir})"
