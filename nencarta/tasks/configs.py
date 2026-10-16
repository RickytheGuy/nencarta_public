import io
from pathlib import Path

import yaml

from nencarta.logger import LOG
from nencarta.workspace import Workspace
from nencarta.core.configs import NencartaConfig
from nencarta.core.defaults import DEFAULT_BATHY_ARGS, DEFAULT_CONFIG, DEFAULT_FLOODMAP_ARGS
from nencarta.core.enumerations import FloodMapMode, Mapper

def _arc_inputs(dem: Path,
                stream_raster: Path,
                land_cover: Path,
                mannings_n_file: Path,
                reanalysis_flow_file: Path,
                river_id: str,
                highflow_field: str,
                baseflow_field: str,
                bathy_args: dict,
                include_baseflow: bool = True,
                ):
    """The input section of an ARC config, with the ARC settings found in ``bathy_args``."""
    params = {
        "#ARC_Inputs": "",
        "DEM_File": dem,
        "Stream_File": stream_raster,
        "LU_Raster_SameRes": land_cover,
        "LU_Manning_n": mannings_n_file,
        "Flow_File": reanalysis_flow_file,
        "Flow_File_ID": river_id,
    }
    if include_baseflow:
        params["Flow_File_BF"] = baseflow_field
    params["Flow_File_QMax"] = highflow_field
    params["Spatial_Units"] = "deg"

    params.update(bathy_args)

    return params

def _fldpln_inputs(workspace: Workspace):
    """The FLDPLN section of a Curve2Flood config."""
    output = {
        "# FLDPLN_Inputs": "",
        "Flow_Direction_File": workspace.flowdir,
        "Filled_DEM_File": workspace.filled_dem,
        "Stream_Info_File": workspace.stream_info_file,
        "FLDPLN_Library": workspace.fldpln_library,
        "reach_id": workspace.configs.stream_id_field,
        "downstream_reach_id": workspace.configs.downstream_id_field,
    }
    output[""] = ""
    return output

def _write_config(config_path: Path, params: dict):
    """
    Write an ARC or Curve2Flood config: YAML for a .yaml path, otherwise tab-separated text. Keys starting with
    '#' are section headers.
    """
    config_path.parent.mkdir(parents=True, exist_ok=True)
    with open(config_path, 'w') as f:
        if config_path.suffix.lower() in {'.yaml', '.yml'}:
            stream = io.StringIO()
            # Convert all Paths in params to strings for YAML serialization
            params_serializable = {key: (str(value) if isinstance(value, (Path, Mapper)) else value) for key, value in params.items()}
            yaml.dump(params_serializable, stream, encoding='utf-8', sort_keys=False, Dumper=yaml.CDumper) # Use CDumper for faster dumping
            # if a line startwith "'#", then strip the quotation and make it a comment
            content = stream.getvalue()
            content = '\n'.join(line if not line.startswith("'#") else f"# {line[2:].split(':')[0][:-1]}" for line in content.split('\n'))
            # Replace lines starting with : or ? with ''
            content = '\n'.join(line if (not line.startswith((':', '?', "'"))) else '' for line in content.split('\n'))
            f.write(content)
        else:
            first_line = True
            for key, value in params.items():
                # Legacy config style: section headers are single tokens on their own line,
                # with a blank line separating sections.
                if isinstance(key, str) and key.startswith("#"):
                    if not first_line:
                        f.write("\n\n")
                    f.write(str(key))
                else:
                    if not first_line:
                        f.write("\n")
                    f.write(f"{key}\t{value}")
                first_line = False

def _has_required_arc_outputs(workspace: Workspace) -> bool:
    """Whether every output ARC is asked for (VDT, curve file, AP database, representative cross sections) exists."""
    configs = workspace.configs
    asked_for = (
        (workspace.VDT_File_Bathy, configs.make_vdt),
        (workspace.Curve_File_Bathy, configs.make_curvefile),
        (workspace.AP_File, configs.make_ap_database),
        (workspace.Representative_Cross_Section_File, configs.make_representative_cross_section_file),
    )
    return all(path.exists() for path, wanted in asked_for if wanted)

def forecast_postfix(configs: NencartaConfig, date: str, hour: str | None) -> str:
    """The end of the flood map names of the forecast issued on ``date`` (and at ``hour``, for NWM)."""
    if configs.streamflow_source.is_nwm():
        return f"_{configs.streamflow_source.file_prefix}_ARC_Flood{configs.floodmap_id}Forecast_{date}_{hour}"
    return f"Forecast_{date}"

def flood_map_paths(workspace: Workspace, postfix: str) -> dict[str, Path]:
    """The flood maps whose names end with ``postfix``, by the mapper config key that writes each."""
    configs = workspace.configs
    suffix = f"{configs.floodmap_id}{postfix}"

    def path(blurb: str, ext: str) -> Path:
        return workspace.file_path(workspace.flood_folder, blurb, ext, 'source', suffix=suffix,
                                   legacy=f"{configs.streamflow_source.file_prefix}_{workspace.FileName}_{blurb}{suffix}.{ext}")

    return {"OutFLD": path('ARC_Flood', 'tif'), "OutSHP": path('ARC_Flood', 'gpkg'), "OutVEL": path('ARC_FloodVEL', 'tif'),
            "OutWSE": path('ARC_FloodWSE', 'tif'), "OutDEP": path('ARC_FloodDepth', 'tif')}
                
def define_configs_for_dem_cleaning(workspace: Workspace) -> Path:
    """
    Write the config of the DEM cleaner's ARC run (initial rating curves) and, with
    use_specified_depth_for_bathy_mask, of its initial flood map at the first specified depth.
    """
    workspace.ARC_Folder.mkdir(parents=True, exist_ok=True)
    workspace.VDT_Folder.mkdir(parents=True, exist_ok=True)
    workspace.flood_folder.mkdir(parents=True, exist_ok=True)
    workspace.dem_updated_folder.mkdir(parents=True, exist_ok=True)

    configs = workspace.configs
    if workspace.ARC_FileName_for_DEM_Cleaner.exists() and not configs.overwrite:
        return workspace.ARC_FileName_for_DEM_Cleaner

    if configs.burn_streams:
        dem = workspace.fixed_dem
    else:
        dem = workspace.assigned_dem

    if configs.move_stream_network_to_thalweg:
        streamlines = workspace.new_StrmShp_matched
        stream_raster = workspace.new_stream_raster
    else:
        streamlines = workspace.DEM_StrmShp
        stream_raster = workspace.STRM_File_Clean

    params = _arc_inputs(
        dem=dem,
        stream_raster=stream_raster,
        land_cover=workspace.LAND_File,
        mannings_n_file=workspace.mannings_n_text_file,
        reanalysis_flow_file=workspace.DEM_Reanalsyis_FlowFile,
        river_id="COMID",
        highflow_field=configs.specified_highflow_field,
        baseflow_field=configs.specified_bathyflow_field,
        bathy_args=configs.bathy_args,
        include_baseflow=True,
    )
    
    params["reach_id"] = configs.stream_id_field
    params["downstream_reach_id"] = configs.downstream_id_field
    if configs.standard_step:
        params["Standard_Step"] = True
    if configs.use_power_laws_for_bathymetry:
        params["drainage_area_field"] = configs.area_km2_field
        params["coefficient_depth"] = configs.coefficient_depth
        params["exponent_depth"] = configs.exponent_depth
        params["coefficient_width"] = configs.coefficient_width
        params["exponent_width"] = configs.exponent_width

    params["#VDT_Output_File_and_CurveFile"] = ""
    params["VDT_Database_NumIterations"] = configs.bathy_args.get("VDT_Database_NumIterations", DEFAULT_BATHY_ARGS["VDT_Database_NumIterations"])
    params["Print_VDT_Database"] = workspace.VDT_File_Initial
    params["Print_Curve_File"] = workspace.Curve_File_Initial
    params["Reach_Average_Curve_File"] = configs.create_reach_average_curve_file

    params["#Mapper Input Data"] = ""
    params["StrmShp_File"] = streamlines
    params["Comid_Flow_File"] = workspace.COMID_Q_File
    params["mapper"] = configs.mapper if not configs.mapper.is_curve2flood_fldpln_mapper() else Mapper.CURVE2FLOOD_KERNEL_WEIGHTED
    params["Bathy_Use_Banks"] = configs.bathy_use_banks
    if configs.flood_waterlc_and_strm_cells:
        params["Flood_WaterLC_and_STRM_Cells"] = True
        params["LAND_WaterValue"] = configs.land_watervalue
    if configs.find_banks_based_on_landcover:
        params["FindBanksBasedOnLandCover"] = True
    if configs.use_specified_depth_for_bathy_mask:
        params["OutFLD"] = workspace.FloodMapFile_Initial
        params["OutSHP"] = workspace.FloodMapFile_Initial_SHP
        params["FloodSpreader_SpecifyDepth"] = configs.specify_depths_for_bathy_mask[0]

    _write_config(workspace.ARC_FileName_for_DEM_Cleaner, params)
    return workspace.ARC_FileName_for_DEM_Cleaner

def define_arc_configs(workspace: Workspace,) -> Path:
    """
    Write the config of ARC's main run (rating curves and bathymetry), which the flood mapper also reads to put the
    bathymetry into the DEM. Without bathymetry, nothing is rewritten if ARC's outputs exist.
    """
    workspace.ARC_Folder.mkdir(parents=True, exist_ok=True)
    workspace.VDT_Folder.mkdir(parents=True, exist_ok=True)
    workspace.flood_folder.mkdir(parents=True, exist_ok=True)

    configs = workspace.configs
    if configs.disable_bathymetry and _has_required_arc_outputs(workspace) and not configs.overwrite:
        LOG.info(f"Domain {workspace.watershed} has required ARC outputs and bathymetry is disabled, skipping ARC.")
        return workspace.ARC_FileName_Bathy

    if configs.move_stream_network_to_thalweg:
        stream_raster = workspace.new_stream_raster
        stream_vector = workspace.new_StrmShp_matched
    else:
        stream_raster = workspace.STRM_File_Clean
        stream_vector = workspace.DEM_StrmShp


    params = _arc_inputs(
        dem=workspace.model_dem,
        stream_raster=stream_raster,
        land_cover=workspace.LAND_File,
        mannings_n_file=workspace.mannings_n_text_file,
        reanalysis_flow_file=workspace.DEM_Reanalsyis_FlowFile,
        river_id="COMID",
        highflow_field=configs.specified_highflow_field,
        baseflow_field=configs.specified_bathyflow_field,
        bathy_args=configs.bathy_args,
        include_baseflow=not configs.disable_bathymetry,
    )
    
    params["# VDT_Output_File_and_CurveFile"] = ""
    params["VDT_Database_NumIterations"] = configs.bathy_args.get("VDT_Database_NumIterations", DEFAULT_BATHY_ARGS["VDT_Database_NumIterations"])
    if configs.make_vdt:
        params["Print_VDT_Database"] = workspace.VDT_File_Bathy
    if configs.make_curvefile:
        params["Print_Curve_File"] = workspace.Curve_File_Bathy
    if configs.make_ap_database:
        params["Print_AP_Database"] = workspace.AP_File
    if configs.make_cross_section_file:
        params["XS_Out_File"] = workspace.Cross_Section_File
    if configs.make_representative_cross_section_file:
        # ARC builds them in this run, from the cross sections it makes the rating curves from
        params["Build_Representative_Cross_Section"] = True
        params["Representative_Cross_Section_File"] = workspace.Representative_Cross_Section_File

    params["Reach_Average_Curve_File"] = configs.create_reach_average_curve_file
    params["StrmShp_File"] = stream_vector
    params["Slope_Low_Percentile"] = configs.slope_low_percentile
    params["Slope_High_Percentile"] = configs.slope_high_percentile
    if configs.standard_step:
        # ARC links the cross sections down the stream network, so it needs each reach's id and the one below it
        params["Standard_Step"] = True
        params["reach_id"] = configs.stream_id_field
        params["downstream_reach_id"] = configs.downstream_id_field

    if not configs.disable_bathymetry:
        params["reach_id"] = configs.stream_id_field
        params["downstream_reach_id"] = configs.downstream_id_field
        if configs.use_power_laws_for_bathymetry:
            params["drainage_area_field"] = configs.area_km2_field
            params["coefficient_depth"] = configs.coefficient_depth or DEFAULT_CONFIG["coefficient_depth"]
            params["exponent_depth"] = configs.exponent_depth or DEFAULT_CONFIG["exponent_depth"]
            params["coefficient_width"] = configs.coefficient_width or DEFAULT_CONFIG["coefficient_width"]
            params["exponent_width"] = configs.exponent_width or DEFAULT_CONFIG["exponent_width"]

        params["# Mapper Input Data"] = ""
        params["Comid_Flow_File"] = workspace.COMID_Q_File
        # Use kernel method for bathymetry creation
        params["mapper"] = configs.mapper if not configs.mapper.is_curve2flood_fldpln_mapper() else Mapper.CURVE2FLOOD_KERNEL_WEIGHTED
        params["TW_MultFact"] = configs.bathy_args.get("TW_MultFact", DEFAULT_FLOODMAP_ARGS["TW_MultFact"])
        params["TopWidthPlausibleLimit"] = configs.bathy_args.get("TopWidthPlausibleLimit", DEFAULT_FLOODMAP_ARGS["TopWidthPlausibleLimit"])
        if not configs.bathy_args.get("Make_Output_GPKG", DEFAULT_FLOODMAP_ARGS["Make_Output_GPKG"]):
            params["Make_Output_GPKG"] = False

        if configs.use_specified_depth_for_bathy_mask:
            params["FloodSpreader_SpecifyDepth"] = configs.specify_depths_for_bathy_mask[-1]
        else:
            if not workspace.bathy_water_mask:
                raise ValueError("Bathy water mask is required when not using specified depth for bathy mask.")
            params["BathyWaterMask"] = workspace.bathy_water_mask

        params["# Bathymetry_Information"] = ""
        params["Bathy_Trap_H"] = configs.bathy_args.get("Bathy_Trap_H", DEFAULT_BATHY_ARGS["Bathy_Trap_H"])
        params["Bathy_Use_Banks"] = configs.bathy_use_banks

        if configs.flood_waterlc_and_strm_cells:
            params["Flood_WaterLC_and_STRM_Cells"] = True
            params["LAND_WaterValue"] = configs.land_watervalue
        if configs.find_banks_based_on_landcover:
            params["FindBanksBasedOnLandCover"] = True

        workspace.bathy_file_folder.mkdir(parents=True, exist_ok=True)
        params["BATHY_Out_File"] = workspace.ARC_BathyFile
        params["AROutBATHY"] = workspace.ARC_BathyFile
        params["FSOutBATHY"] = workspace.FS_BathyFile
    
    LOG.info(f"Writing ARC config to {workspace.ARC_FileName_Bathy}.")
    _write_config(workspace.ARC_FileName_Bathy, params)
    
    return workspace.ARC_FileName_Bathy

def define_mapper_configs(workspace: Workspace, flow_file: Path) -> Path:
    """
    Write the flood mapper's config for one flow file: a forecast, a user flow file or a return period. Its flood
    maps are named after the flows.
    """
    workspace.flood_folder.mkdir(parents=True, exist_ok=True)
    
    configs = workspace.configs
    config_ext = 'yaml' if configs.use_yaml else 'txt'
    params = {'#ARC_Inputs': ''}
    params["DEM_File"] = workspace.model_dem if configs.disable_bathymetry else workspace.FS_BathyFile
    params["Stream_File"] = workspace.new_stream_raster if configs.move_stream_network_to_thalweg else workspace.STRM_File_Clean
    params["LU_Manning_n"] = workspace.mannings_n_text_file

    params["# VDT_Output_File_and_CurveFile"] = ""

    # The mapper reads the VDT when it is given, and otherwise the curve file
    if configs.make_vdt:
        params["Print_VDT_Database"] = workspace.VDT_File_Bathy
    if configs.make_curvefile:
        params["Print_Curve_File"] = workspace.Curve_File_Bathy
    
    params["# Mapper Input Data"] = ""
    params["StrmShp_File"] = workspace.new_StrmShp_matched if configs.move_stream_network_to_thalweg else workspace.DEM_StrmShp
    params["Comid_Flow_File"] = flow_file
    params["mapper"] = configs.mapper
    if configs.mapper.is_curve2flood_fldpln_mapper():
        params.update(_fldpln_inputs(workspace))
    params["TW_MultFact"] = configs.floodmap_args.get("TW_MultFact", DEFAULT_FLOODMAP_ARGS["TW_MultFact"])
    params["TopWidthPlausibleLimit"] = configs.floodmap_args.get("TopWidthPlausibleLimit", DEFAULT_FLOODMAP_ARGS["TopWidthPlausibleLimit"])

    if configs.flood_waterlc_and_strm_cells or configs.make_velocity_maps:
        params["Flood_WaterLC_and_STRM_Cells"] = configs.flood_waterlc_and_strm_cells
        params["LU_Raster_SameRes"] = workspace.LAND_File
        params["LAND_WaterValue"] = configs.land_watervalue
    if configs.flood_bathy_mask:
        params["Flood_Bathy_Mask"] = True
        params["BathyWaterMask"] = workspace.bathy_water_mask

    if configs.floodmap_mode == FloodMapMode.FORECAST:
        if configs.forensic_forecast_date:
            postfix = forecast_postfix(configs, configs.forensic_forecast_date, configs.forensic_forecast_hour)
        else:
            postfix = forecast_postfix(configs, configs.forecastdate, configs.forecasthour)
        config_path = workspace.file_path(workspace.ARC_Folder, 'ARC_Input_FloodForecast', config_ext, 'source',
                                          legacy=f"{configs.streamflow_source.file_prefix}_ARC_Input_{workspace.FileName}_FloodForecast.{config_ext}")
    elif configs.floodmap_mode == FloodMapMode.USER:
        postfix = f"_{flow_file.stem}"
        if not flow_file.exists():
            LOG.error(f"User provided flow file does not exist: {flow_file}")
            raise FileNotFoundError(f"User provided flow file does not exist: {flow_file}")
        config_path = workspace.file_path(workspace.ARC_Folder, 'ARC_Input', config_ext, 'source', suffix=f"{configs.floodmap_id}{postfix}",
                                          legacy=f"{configs.streamflow_source.file_prefix}_ARC_Input_{workspace.FileName}{configs.floodmap_id}{postfix}.{config_ext}")
    elif configs.floodmap_mode == FloodMapMode.RETURN_PERIOD:
        return_period = next((rp for rp in configs.return_periods if workspace.return_period_flow_file(rp) == Path(flow_file)), None)
        postfix = f"_rp{return_period}" if return_period is not None else f"_{flow_file.stem.rsplit('_')[-1]}"
        config_path = workspace.file_path(workspace.ARC_Folder, 'ARC_Input', config_ext, 'source', suffix=postfix,
                                          legacy=f"{configs.streamflow_source.file_prefix}_ARC_Input_{workspace.FileName}{postfix}.{config_ext}")
    else:
        LOG.error(f"Invalid floodmap_mode: {configs.floodmap_mode}")
        raise ValueError(f"Invalid floodmap_mode: {configs.floodmap_mode}")

    flood_maps = flood_map_paths(workspace, postfix)
    params["OutFLD"] = flood_maps["OutFLD"]

    if configs.floodmap_args.get('Make_Output_GPKG', DEFAULT_FLOODMAP_ARGS["Make_Output_GPKG"]):
        params["OutSHP"] = flood_maps["OutSHP"]
    else:
        params["Make_Output_GPKG"] = False

    if configs.make_velocity_maps:
        params["OutVEL"] = flood_maps["OutVEL"]

    if configs.make_wse_maps:
        params["OutWSE"] = flood_maps["OutWSE"]

    if configs.make_depth_maps:
        params["OutDEP"] = flood_maps["OutDEP"]

    LOG.info(f"Writing Mapper config to {config_path}.")
    _write_config(config_path, params)
    return config_path
