from pathlib import Path

import yaml

from nencarta.api.enumerations import FloodMapMode, Mapper, StreamflowSource
from nencarta.logger import LOG
from nencarta.workspace import Workspace
from nencarta.api.configs import NencartaConfig

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
    # Keep key names + ordering aligned with the legacy configs produced by main.py.
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

    optional_bathy_keys = [
        "X_Section_Dist",
        "Degree_Manip",
        "Degree_Interval",
        "Low_Spot_Range",
        "Str_Limit_Val",
        "Gen_Dir_Dist",
        "Gen_Slope_Dist",
        "Stream_Slope_Method",
    ]

    params.update({
        key: bathy_args[key]
        for key in optional_bathy_keys
        if key in bathy_args
    })

    return params

def _fldpln_inputs(workspace: Workspace):
    return {
        "# FLDPLN_Inputs": "",
        "Flow_Direction_File": workspace.flowdir,
        "Filled_DEM_File": workspace.filled_dem,
        "Stream_Info_File": workspace.stream_info_file,
        "FLDPLN_Library": workspace.fldpln_library,
        "FSOutBATHY": workspace.FS_BathyFile,
    }

def _write_config(config_path: Path, params: dict):
    config_path.parent.mkdir(parents=True, exist_ok=True)
    with open(config_path, 'w') as f:
        if config_path.suffix.lower() in {'.yaml', '.yml'}:
            # Convert all Paths in params to strings for YAML serialization
            params_serializable = {key: (str(value) if isinstance(value, (Path, Mapper)) else value) for key, value in params.items()}
            yaml.dump(params_serializable, f, encoding='utf-8', sort_keys=False)
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
    required_files = [workspace.VDT_File_Bathy]
    if workspace.Curve_File_Bathy:
        required_files.append(workspace.Curve_File_Bathy)
    if workspace.AP_File:
        required_files.append(workspace.AP_File)
    return all(path.exists() for path in required_files)
                
def define_configs_for_dem_cleaning(workspace: Workspace) -> Path:
    workspace.ARC_Folder.mkdir(parents=True, exist_ok=True)
    workspace.VDT_Folder.mkdir(parents=True, exist_ok=True)
    workspace.flood_folder.mkdir(parents=True, exist_ok=True)
    workspace.dem_updated_folder.mkdir(parents=True, exist_ok=True)

    configs = workspace.configs
    if workspace.ARC_FileName_for_DEM_Cleaner.exists() and not configs.process_stream_network:
        return workspace.ARC_FileName_for_DEM_Cleaner

    params = _arc_inputs(
        dem=workspace.assigned_dem,
        stream_raster=workspace.STRM_File_Clean,
        land_cover=workspace.LAND_File,
        mannings_n_file=workspace.mannings_n_text_file,
        reanalysis_flow_file=workspace.DEM_Reanalsyis_FlowFile,
        river_id="COMID",
        highflow_field=configs.specified_highflow_field,
        baseflow_field=configs.specified_bathyflow_field,
        bathy_args=configs.bathy_args,
        include_baseflow=True,
    )

    params["#VDT_Output_File_and_CurveFile"] = ""
    params["VDT_Database_NumIterations"] = 30
    params["Print_VDT_Database"] = workspace.VDT_File_Initial
    params["Print_Curve_File"] = workspace.Curve_File_Initial
    params["Reach_Average_Curve_File"] = configs.create_reach_average_curve_file

    params["#Mapper Input Data"] = ""
    params["StrmShp_File"] = workspace.DEM_StrmShp
    params["Comid_Flow_File"] = workspace.COMID_Q_File
    params["mapper"] = configs.mapper
    params["FS_ADJUST_FLOW_BY_FRACTION"] = 1.0
    params["Bathy_Use_Banks"] = configs.bathy_use_banks
    if configs.flood_waterlc_and_strm_cells:
        params["LAND_WaterValue"] = configs.land_watervalue
    if configs.find_banks_based_on_landcover:
        params["FindBanksBasedOnLandCover"] = True
    if configs.use_specified_depth_for_bathy_mask:
        if configs.mapper.is_curve2flood_fldpln_mapper():
            params.update(_fldpln_inputs(workspace))
        params["OutFLD"] = workspace.FloodMapFile_Initial
        params["OutSHP"] = workspace.FloodMapFile_Initial_SHP
        params["FloodSpreader_SpecifyDepth"] = configs.specify_depths_for_bathy_mask[0]

    _write_config(workspace.ARC_FileName_for_DEM_Cleaner, params)
    return workspace.ARC_FileName_for_DEM_Cleaner

def define_arc_configs(workspace: Workspace,) -> Path:
    workspace.ARC_Folder.mkdir(parents=True, exist_ok=True)
    workspace.VDT_Folder.mkdir(parents=True, exist_ok=True)
    workspace.flood_folder.mkdir(parents=True, exist_ok=True)

    configs = workspace.configs
    if configs.disable_bathymetry and _has_required_arc_outputs(workspace) and not configs.process_stream_network:
        LOG.info(f"Domain {workspace.watershed} has required ARC outputs and bathymetry is disabled, skipping ARC.")
        return workspace.ARC_FileName_Bathy

    params = _arc_inputs(
        dem=workspace.assigned_dem,
        stream_raster=workspace.STRM_File_Clean,
        land_cover=workspace.LAND_File,
        mannings_n_file=workspace.mannings_n_text_file,
        reanalysis_flow_file=workspace.DEM_Reanalsyis_FlowFile,
        river_id="COMID",
        highflow_field=configs.specified_highflow_field,
        baseflow_field=configs.specified_bathyflow_field,
        bathy_args=configs.bathy_args,
        include_baseflow=not configs.disable_bathymetry,
    )
    
    params["#VDT_Output_File_and_CurveFile"] = ""
    params["VDT_Database_NumIterations"] = configs.bathy_args.get("VDT_Database_NumIterations", 30)
    params["Print_VDT_Database"] = workspace.VDT_File_Bathy
    if configs.make_curvefile:
        params["Print_Curve_File"] = workspace.Curve_File_Bathy
    if configs.make_ap_database:
        params["Print_AP_Database"] = workspace.AP_File
        # TODO
    # params["XS_Out_File"] = workspace.VDT_Folder / f"{workspace.FileName}_XS.csv"
    params["Reach_Average_Curve_File"] = configs.create_reach_average_curve_file

    if not configs.disable_bathymetry:
        params["#Mapper Input Data"] = ""
        params["Comid_Flow_File"] = workspace.COMID_Q_File
        params["StrmShp_File"] = workspace.DEM_StrmShp
        params["mapper"] = configs.mapper
        params["FS_ADJUST_FLOW_BY_FRACTION"] = configs.bathy_args.get("FS_ADJUST_FLOW_BY_FRACTION", 1.0)
        params["OutFLD"] = workspace.FloodMapFile_Bathy
        params["OutSHP"] = workspace.FloodMapFile_Bathy_SHP
        params["TopWidthDistanceFactor"] = configs.bathy_args.get("TopWidthDistanceFactor", 1.5)
        params["TW_MultFact"] = configs.bathy_args.get("TW_MultFact", 1.5)
        params["TopWidthPlausibleLimit"] = configs.bathy_args.get("TopWidthPlausibleLimit", 2000)
        if not configs.bathy_args.get("Make_Output_GPKG", True):
            params["Make_Output_GPKG"] = False
        if configs.mapper.is_curve2flood_fldpln_mapper():
            params.update(_fldpln_inputs(workspace))

        if configs.use_specified_depth_for_bathy_mask:
            params["FloodSpreader_SpecifyDepth"] = configs.specify_depths_for_bathy_mask[-1]
        else:
            if not workspace.bathy_water_mask:
                raise ValueError("Bathy water mask is required when not using specified depth for bathy mask.")
            params["BathyWaterMask"] = workspace.bathy_water_mask

        params["#Bathymetry_Information"] = ""
        params["Bathy_Trap_H"] = configs.bathy_args.get("Bathy_Trap_H", 0.2)
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
    workspace.flood_folder.mkdir(parents=True, exist_ok=True)
    
    configs = workspace.configs
    params = {'#ARC_Inputs': ''}
    if configs.disable_bathymetry:
        params["DEM_File"] = workspace.assigned_dem
    else:
        params["DEM_File"] = workspace.FS_BathyFile
    params["Stream_File"] = workspace.STRM_File_Clean
    params["LU_Manning_n"] = workspace.mannings_n_text_file

    params["#VDT_Output_File_and_CurveFile"] = ""
    
    params["Print_VDT_Database"] = workspace.VDT_File_Bathy
    if configs.make_curvefile:
        params["Print_Curve_File"] = workspace.Curve_File_Bathy
    
    params["#Mapper Input Data"] = ""
    params["StrmShp_File"] = workspace.DEM_StrmShp
    params["Comid_Flow_File"] = flow_file
    params["mapper"] = configs.mapper
    if configs.mapper.is_curve2flood_fldpln_mapper():
        params.update(_fldpln_inputs(workspace))
    params["FS_ADJUST_FLOW_BY_FRACTION"] = configs.floodmap_args.get("FS_ADJUST_FLOW_BY_FRACTION", 1.0)
    params["TW_MultFact"] = configs.floodmap_args.get("TW_MultFact", 1.5)
    params["TopWidthPlausibleLimit"] = configs.floodmap_args.get("TopWidthPlausibleLimit", 6000)

    if configs.flood_waterlc_and_strm_cells or configs.make_velocity_maps or True: # TODO
        params["Flood_WaterLC_and_STRM_Cells"] = configs.flood_waterlc_and_strm_cells
        params["LU_Raster_SameRes"] = workspace.LAND_File
        params["LAND_WaterValue"] = configs.land_watervalue


    if configs.floodmap_mode == FloodMapMode.FORECAST:
        if configs.streamflow_source.is_nwm():
            # create the end of the file name that describes the forecast
            if configs.forensic_forecast_date != None:
                ending_of_forecast_file = f'Forecast_{configs.forensic_forecast_date}_{configs.forensic_forecast_hour}'
            elif configs.forecastdate != None:
                ending_of_forecast_file = f'Forecast_{configs.forecastdate}_{configs.forecasthour}'
            # rename the forecast of the extent raster based upon the type of NWM forecast we are using
            postfix = f"_{configs.streamflow_source}_ARC_Flood{configs.floodmap_id}{ending_of_forecast_file}"
        elif configs.streamflow_source == StreamflowSource.GEOGLOWS:
            # create the end of the file name that describes the forecast
            if configs.forensic_forecast_date != None:
                ending_of_forecast_file = f'Forecast_{configs.forensic_forecast_date}'
            elif configs.forecastdate != None:
                ending_of_forecast_file = f'Forecast_{configs.forecastdate}'
            postfix = ending_of_forecast_file
            config_path = workspace.ARC_Folder / f"{configs.streamflow_source}_ARC_Input_{workspace.FileName}_FloodForecast.{'yaml' if configs.use_yaml else 'txt'}"
    elif configs.floodmap_mode == FloodMapMode.USER:
        postfix = f"_{flow_file.stem}"
        if not flow_file.exists():
            LOG.error(f"User provided flow file does not exist: {flow_file}")
            raise FileNotFoundError(f"User provided flow file does not exist: {flow_file}")
        config_path = workspace.ARC_Folder / f"{configs.streamflow_source}_ARC_Input_{workspace.FileName}{configs.floodmap_id}{postfix}.{'yaml' if configs.use_yaml else 'txt'}"
    elif configs.floodmap_mode == FloodMapMode.RETURN_PERIOD:
        postfix = f"_{flow_file.stem.rsplit('_')[-1]}"
        config_path = workspace.ARC_Folder / f"{configs.streamflow_source}_ARC_Input_{workspace.FileName}{postfix}.{'yaml' if configs.use_yaml else 'txt'}"
    else:
        LOG.error(f"Invalid floodmap_mode: {configs.floodmap_mode}")
        raise ValueError(f"Invalid floodmap_mode: {configs.floodmap_mode}")

    floodmap_path = workspace.flood_folder / f"{configs.streamflow_source}_{workspace.FileName}_ARC_Flood{configs.floodmap_id}{postfix}.tif"
    params["OutFLD"] = floodmap_path

    if configs.floodmap_args.get('Make_Output_GPKG', True):
        floodmap_vector = workspace.flood_folder / f"{configs.streamflow_source}_{workspace.FileName}_ARC_Flood{configs.floodmap_id}{postfix}.gpkg"
        params["OutSHP"] = floodmap_vector
    else:
        params["Make_Output_GPKG"] = False

    if configs.make_velocity_maps:
        velocity_path = workspace.flood_folder / f"{configs.streamflow_source}_{workspace.FileName}_ARC_FloodVEL{configs.floodmap_id}{postfix}.tif"
        params["OutVEL"] = velocity_path

    if configs.make_wse_maps:
        wse_path = workspace.flood_folder / f"{configs.streamflow_source}_{workspace.FileName}_ARC_FloodWSE{configs.floodmap_id}{postfix}.tif"
        params["OutWSE"] = wse_path

    if configs.make_depth_maps:
        depth_path = workspace.flood_folder / f"{configs.streamflow_source}_{workspace.FileName}_ARC_FloodDepth{configs.floodmap_id}{postfix}.tif"
        params["OutDEP"] = depth_path

    LOG.info(f"Writing Mapper config to {config_path}.")
    _write_config(config_path, params)
    return config_path