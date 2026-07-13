import re
from pathlib import Path
from functools import cache
from datetime import datetime, timedelta

import pandas as pd

from nencarta.logger import LOG
from nencarta.core.vector import Vector
from nencarta.core.configs import NencartaConfig
from nencarta.workspace import Workspace
from nencarta.core.floodmapper_output import FloodMapperOutput
from nencarta import Download_Process_ForecastData as ForecastFlows

def make_flow_file_from_forecast(workspace: Workspace) -> list[Path]:
    # now lets download the forecast streamflows
    #Forecast flow data from GeoGLOWS
    # parquet_file_from_geoglows = 'v2-model-table.parquet'     #http://geoglows-v2.s3-website-us-west-2.amazonaws.com/#tables/
    configs = workspace.configs
    rivids = Vector(workspace.DEM_StrmShp, not workspace.configs.parallel).to_geopandas(columns=[configs.stream_id_field])[configs.stream_id_field].astype(int).tolist()
    
    forecastdate = configs.forensic_forecast_date
    forecasthour = configs.forensic_forecast_hour
    if forecastdate:
        LOG.info(f"Using forensic forecast date: {forecastdate}")
        if configs.streamflow_source.is_nwm():
            LOG.info(f"Using forensic forecast hour: {forecasthour}")
        else:
            forecasthour = None
        
        if configs.streamflow_source.is_nwm():
            flow_file_name = f'{workspace.FileName}_{str(forecastdate)}_{forecasthour}_{configs.streamflow_source}_forecast.csv'
        else:
            flow_file_name = f'{workspace.FileName}_{str(forecastdate)}_{configs.streamflow_source}_forecast.csv'
            
        try:
            ForecastFlowFile = workspace.FLOW_Folder / flow_file_name
            if not ForecastFlowFile.exists() or configs.process_stream_network:
                ForecastFlows.Process_and_Write_Forecast_Data(forecastdate, forecasthour, rivids, ForecastFlowFile, configs.streamflow_source, configs.nwm_api_key)
        except Exception as e:
            LOG.error('Could not process forensic forecast streamflow download, please check your date or try again later...')
            raise e
    else:
        # cycle through today through 12 days ago to find the most recent day with a forecast
        found = False
        for fd in range(0,13):
            for fh in range(0,24):
                try:
                    ForecastFlowFile = workspace.FLOW_Folder / flow_file_name
                    forecastdate, forecasthour = ForecastFlows.Get_Date_For_Forecast(fd, fh, configs.streamflow_source) 
                    LOG.info(f"Attempting to download forecast for date: {forecastdate} and hour: {forecasthour}...")
                    # we only need the forecast date for GEOGLOWS, for NWM we need the forecast hour as well               
                    if configs.streamflow_source.is_nwm():
                        flow_file_name = f'{workspace.FileName}_{str(forecastdate)}_{forecasthour}_{configs.streamflow_source}_forecast.csv'
                    else:
                        flow_file_name = f'{workspace.FileName}_{str(forecastdate)}_{configs.streamflow_source}_forecast.csv'

                    ForecastFlowFile = workspace.FLOW_Folder / flow_file_name
                    if not ForecastFlowFile.exists() or configs.process_stream_network:
                        ForecastFlows.Process_and_Write_Forecast_Data(forecastdate, forecasthour, rivids, ForecastFlowFile, configs.streamflow_source, configs.nwm_api_key)

                    if configs:
                        configs.forecastdate = forecastdate
                        configs.forecasthour = forecasthour

                    found = True
                    break
                except Exception as e:
                    LOG.error(f'Could not process forecast, moving back another day.. ({e})')
            if found:
                break  # break outer

    LOG.info(f'Forecast data save here: {ForecastFlowFile}')    
    return [ForecastFlowFile]

def assign_user_flow_files(configs: NencartaConfig,) -> list[Path]:
    flow_files = []
    if configs.user_flow_files:
        for file in configs.user_flow_files:
            flow_files.append(Path(file))
    return flow_files

def remove_old_forecast_files(workspace: Workspace, forecast_floodmap_basename: str):
    # check and see if forecasts past a specified date exist and if so, delete them
    for filename in workspace.flood_folder.iterdir():
        if not filename.name.startswith(forecast_floodmap_basename[:-12]) or not filename.is_file():
            continue

        # Regular expression to extract the date
        date_pattern = re.compile(r'\d{8}')
        match = date_pattern.search(filename.name)

        if not match:
            LOG.warning("No valid date found in the filename.")
            continue

        # Extracted date string
        date_str = match.group()
        
        # Convert to datetime object
        file_date = datetime.strptime(date_str, '%Y%m%d')
        
        # Calculate the date
        date_threshold = datetime.now() - timedelta(days=workspace.configs.age_of_forecast_days)
        
        # Check if the file is older than the specified number of days
        if file_date <= date_threshold:
            # If the file is as old or older than the specified number of days, delete the file
            filename.unlink()
            LOG.info(f"File {filename} has been deleted.")
        else:
            LOG.info(f"File {filename} is not old enough to be deleted.")

def make_return_period_flow_file(workspace: Workspace) -> list[Path]:
    # make a flow file with return period flows for each reach based on the specified field in the source flow file
    configs = workspace.configs
    if configs.streamflow_source.is_nwm():
        LOG.error("Return period flow file generation is not currently supported for NWM streamflow source.")
        raise NotImplementedError("Return period flow file generation is not currently supported for NWM streamflow source.")
    
    files = []
    for rp in configs.return_periods:
        return_period_flow_file = workspace.FLOW_Folder / f'{workspace.FileName}_rp{rp}.csv'
        files.append(return_period_flow_file)

        if return_period_flow_file.exists() and not configs.process_stream_network:
            LOG.info(f"{return_period_flow_file} already exists and we aren't making it again...")
            continue
        
        LOG.info(f"Creating return period flow file for rp{rp} at {return_period_flow_file}...")
        pd.read_csv(workspace.DEM_Reanalsyis_FlowFile)[['COMID', f"rp{rp}"]].to_csv(return_period_flow_file, index=False)

    return files