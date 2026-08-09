from pathlib import Path

import pandas as pd

from nencarta.core.model_config import ModelConfig
from nencarta.logger import LOG
from nencarta.core.configs import NencartaConfig
from nencarta.workspace import Workspace
from nencarta.core.enumerations import FloodMapMode, StreamflowSource

def get_fist_inputs(workspace: Workspace, flow_file: Path) -> tuple:
    OutProjection = "EPSG:4269"
    configs = workspace.configs
    tasks = []
    # SEED file for creating a GEOJSON for FIST
    if configs.floodmap_mode == FloodMapMode.FORECAST:
        # There will be only one file being used here for all forecasts
        SEED_File = workspace.FIST_Folder / f'{workspace.FileName}_Seed.parquet'
    elif configs.floodmap_mode == FloodMapMode.USER:
        SEED_File = workspace.FIST_Folder / f"{workspace.FileName}_{flow_file.name.rsplit('.', 1)[0]}_Seed.parquet"

    if SEED_File.exists() and not configs.overwrite:
        LOG.info(f"Seed file {SEED_File} already exists for domain {workspace.watershed}, skipping FIST GeoJSON generation.")
        return tasks


    streamflow_forecast_df = pd.read_csv(flow_file)
    streamflow_columns = streamflow_forecast_df.select_dtypes(include=['float']).columns.tolist()

    # grab the stream id column, which should be the first column in the flow file and should be named 'rivid' or 'comid' depending on the source of streamflow data
    id_column_name = streamflow_forecast_df.columns[0]

    for streamflow_column in streamflow_columns:            
        streamflow_forecast_filtered_df = streamflow_forecast_df[[id_column_name, streamflow_column]]
        if configs.floodmap_mode == FloodMapMode.FORECAST:
            # if it is a forecast, decifer the name of the forecast based upon the type of streamflow data and the presence of a forensic forecast date, and then name the geojson accordingly
            if configs.forensic_forecast_date is not None and configs.streamflow_source == StreamflowSource.GEOGLOWS:
                GeoJSON_File = workspace.FIST_Folder / f"{workspace.FileName}_{configs.forensic_forecast_date}_{streamflow_column}.geojson"
            elif configs.forensic_forecast_date is None and configs.streamflow_source == StreamflowSource.GEOGLOWS:
                GeoJSON_File = workspace.FIST_Folder / f"{workspace.FileName}_{configs.forecastdate}_{streamflow_column}.geojson"
            elif configs.forensic_forecast_date is not None and configs.streamflow_source.is_nwm():
                GeoJSON_File = workspace.FIST_Folder / f"{workspace.FileName}_{configs.forensic_forecast_date}_{configs.forensic_forecast_hour}_{streamflow_column}.geojson"
            elif configs.forensic_forecast_date is None and configs.streamflow_source.is_nwm():
                GeoJSON_File = workspace.FIST_Folder / f"{workspace.FileName}_{configs.forecastdate}_{configs.forecasthour}_{streamflow_column}.geojson"

        elif configs.floodmap_mode == FloodMapMode.USER:
            GeoJSON_File = workspace.FIST_Folder / f"{workspace.FileName}_{flow_file.name.rsplit('.', 1)[0]}_{streamflow_column}.geojson"

        # Always regenerate the FIST GeoJSON so reruns pick up CRS and
        # geometry fixes instead of silently reusing a stale output file.
        if GeoJSON_File.exists():
            GeoJSON_File.unlink()
        
        tasks.append((
            workspace.VDT_File_Bathy,
            workspace.STRM_File_Clean,
            GeoJSON_File,
            OutProjection,
            workspace.DEM_StrmShp,
            configs.streamflow_source.upstream_id,
            configs.streamflow_source.downstream_id,
            SEED_File,
            True, 
            None,
            streamflow_forecast_filtered_df
        ))

    return tasks

def run_fist(model_config: ModelConfig, flag = None) -> None:
    from arc import Run_Main_Curve_to_GEOJSON_Program_Stream_Vector
    for args in model_config.fist_inputs:
        Run_Main_Curve_to_GEOJSON_Program_Stream_Vector(*args)
