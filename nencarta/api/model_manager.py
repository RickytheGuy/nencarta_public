import os
import json
import subprocess
from contextlib import nullcontext
from functools import partial
import multiprocessing as mp
from pathlib import Path

try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

import tqdm
import pandas as pd
from osgeo import gdal

gdal.UseExceptions()

from nencarta.api.enumerations import Mapper, FloodMapMode, StreamflowSource
from nencarta.api.domain import Domain
from nencarta.processing_tracker import Timer
from nencarta.logger import LOG
from nencarta._constants import FLOODSPREADER_PATH

def _run_arc(config: str, **kwargs):
    from arc import Arc
    Arc(config, **kwargs).run()

def _run_mapper(config: str, mapper: Mapper, **kwargs):
    if mapper == Mapper.FLOODSPREADER:
        if not os.path.exists(FLOODSPREADER_PATH):
            raise FileNotFoundError(f"FloodSpreader script not found at {FLOODSPREADER_PATH}. Please ensure it is included in the nencarta package.")
        call_mapper = f'python "{FLOODSPREADER_PATH}" {config}'
        subprocess.call(call_mapper, shell=True)
    else:
        from curve2flood import Curve2Flood_MainFunction
        Curve2Flood_MainFunction(config, **kwargs)

def _run_fist(args):
    from arc import Run_Main_VDT_to_GEOJSON_Program_Stream_Vector
    Run_Main_VDT_to_GEOJSON_Program_Stream_Vector(*args)

def has_required_arc_outputs(domain: Domain) -> bool:
    required_files = [domain.vdt]
    if domain.curve_file:
        required_files.append(domain.curve_file)
    if domain.ap_file:
        required_files.append(domain.ap_file)
    return all(os.path.exists(path) for path in required_files)

def Create_Go_Consequence_GeoJSON(Consequences_JSON_Path, Forecast_Flood_Depth_Raster_Name, Consequences_Output_GPKG_File):
    """
    Write a go-consequences config JSON using NSIAPI structures,
    one 'depth' hazard GeoTIFF, and a GPKG results writer.

    Args:
        Consequences_GeoJSON_Path (str): Full path (including filename) to write the JSON config on your host.
        FloodDepthFile (str): Container-visible path to the flood-depth GeoTIFF (e.g., '/data/.../file.tif').
        Consequences_GeoJSON_File (str): Container-visible output GPKG path (e.g., '/data/.../file.gpkg').

    Returns:
        None
    """

    config = {
        "structure_provider_info": {
            "structure_provider_type": "NSIAPI"
        },
        "hazard_provider_info": {
            "hazards": [
                {
                    "hazard_parameter_type": "depth",
                    "hazard_provider_file_path": f"/data/FloodMap/{Forecast_Flood_Depth_Raster_Name}"
                }
            ]
        },
        "results_writer_info": {
            "results_writer_type": "GPKG",
            "output_file_path": f"/data/Consequences/{Consequences_Output_GPKG_File}"
        }
    }

    with open(Consequences_JSON_Path, "w", encoding="utf-8") as f:
        json.dump(config, f, indent=4)

    return


class ModelManager:
    def __init__(self, 
                 domains: Domain | list[Domain], 
                 timer: Timer = None,             
                 processes: int = 1,):
        self.domains = domains if isinstance(domains, list) else [domains]
        # Sort domains by their DEM size, largest to smallest, to optimize memory usage
        self.domains.sort(key=lambda d: d.get_priority())

        if timer is None:
            timer = Timer()
        self.timer = timer
        if processes < 0:
            processes = mp.cpu_count() + processes

        if processes < 1:
            raise ValueError("Number of processes must be at least 1.")
        
        self.processes = processes
        self.pool = None

    def __enter__(self) -> Self:
        if self.processes > 1:
            self.pool = mp.Pool(self.processes)
        return self
    
    def __exit__(self, exc_type, exc, tb):
        if self.pool:
            self.pool.close()
            self.pool.join()

    def iter(self, func, items, chunksize=None):
        if self.pool:
            if chunksize is None:
                chunksize = max(1, len(items) // (self.processes * 4))
            return self.pool.imap_unordered(func, items, chunksize=chunksize)
        else:
            for item in items:
                yield func(item)

    def get_arc_tasks(self, overwrite=False) -> list[str]:
        tasks = []
        for domain in self.domains:
            if domain.arc_config and (not os.path.exists(domain.vdt) or overwrite):
                tasks.append(domain.arc_config)
        return tasks
    
    def get_bathy_arc_tasks(self, 
                        bathymetry_disabled: bool,
                        overwrite: bool) -> list[str]:
        arc_tasks = []

        for domain in self.domains:
            if not domain.arc_config:
                continue

            if bathymetry_disabled and has_required_arc_outputs(domain) and not overwrite:
                LOG.info(f"Domain {domain.name} has required ARC outputs and bathymetry is disabled, skipping ARC.")
                continue

            LOG.info(f"Creating {domain.arc_bathy} for domain {domain.name}.")
            arc_tasks.append(domain.arc_config)

        return arc_tasks

    def get_bathy_mapper_tasks(self, 
                        bathymetry_disabled: bool,
                        overwrite: bool) -> list[str]:
        mapper_tasks = []
        for domain in self.domains:
            if not domain.arc_config or bathymetry_disabled:
                continue

            if not os.path.exists(domain.dem_with_bathymetry) or overwrite:
                LOG.info(f"Creating {domain.dem_with_bathymetry} for domain {domain.name}.")
                mapper_tasks.append(domain.arc_config)

        return mapper_tasks

    def get_floodmap_tasks(self, overwrite=False) -> list[str]:
        tasks = []
        for domain in self.domains:
            if domain.floodmap_configs:
                for index, config in enumerate(domain.floodmap_configs):
                    if config and (not os.path.exists(domain.flood_maps[index]) or overwrite):
                        tasks.append(config)
        return tasks
    
    def get_FIST_tasks(self, 
                       floodmap_mode: FloodMapMode,
                       streamflow_source: StreamflowSource,
                       stream_id_field: str,
                       downstream_id_field: str,
                       forensic_forecast_date: str | None,
                       forensic_forecast_hour: str | None,
                       forecastdate: str | None,
                       forecasthour: str | None,
                       overwrite=False) -> list[str]:
        OutProjection = "EPSG:4269"
        tasks = []
        for domain in self.domains:
            flow_files = domain.flood_flow_files
            for flow_file in flow_files:
                # SEED file for creating a GEOJSON for FIST
                if floodmap_mode == FloodMapMode.FORECAST:
                    # There will be only one forecast file being used here
                    SEED_File = os.path.join(domain.folder_structure.FIST_Folder, f'{domain.dem_basename}_Seed.shp') 
                elif floodmap_mode == FloodMapMode.USER:
                    SEED_File = os.path.join(domain.folder_structure.FIST_Folder, f"{domain.dem_basename}_{os.path.basename(flow_file).rsplit('.', 1)[0]}_Seed.shp")

                if os.path.exists(SEED_File) and not overwrite:
                    LOG.info(f"Seed file {SEED_File} already exists for domain {domain.name}, skipping FIST GeoJSON generation.")
                    continue


                streamflow_forecast_df = pd.read_csv(flow_file)
                streamflow_columns = streamflow_forecast_df.select_dtypes(include=['float']).columns.tolist()

                # grab the stream id column, which should be the first column in the flow file and should be named 'rivid' or 'comid' depending on the source of streamflow data
                id_column_name = streamflow_forecast_df.columns[0]

                for streamflow_column in streamflow_columns:            
                    streamflow_forecast_filtered_df = streamflow_forecast_df[[id_column_name, streamflow_column]]
                    if floodmap_mode == FloodMapMode.FORECAST:
                        # if it is a forecast, decifer the name of the forecast based upon the type of streamflow data and the presence of a forensic forecast date, and then name the geojson accordingly
                        if forensic_forecast_date is not None and streamflow_source == StreamflowSource.GEOGLOWS:
                            GeoJSON_File = os.path.join(domain.folder_structure.FIST_Folder, f"{domain.dem_basename}_{forensic_forecast_date}_{streamflow_column}.geojson") 
                        elif forensic_forecast_date is None and streamflow_source == StreamflowSource.GEOGLOWS:
                            GeoJSON_File = os.path.join(domain.folder_structure.FIST_Folder, f"{domain.dem_basename}_{forecastdate}_{streamflow_column}.geojson") 
                        elif forensic_forecast_date is not None and streamflow_source.is_nwm():
                            GeoJSON_File = os.path.join(domain.folder_structure.FIST_Folder, f"{domain.dem_basename}_{forensic_forecast_date}_{forensic_forecast_hour}_{streamflow_column}.geojson")
                        elif forensic_forecast_date is None and streamflow_source.is_nwm():
                            GeoJSON_File = os.path.join(domain.folder_structure.FIST_Folder, f"{domain.dem_basename}_{forecastdate}_{forecasthour}_{streamflow_column}.geojson")

                    elif floodmap_mode == FloodMapMode.USER:
                        GeoJSON_File = os.path.join(domain.folder_structure.FIST_Folder, f"{domain.dem_basename}_{os.path.basename(flow_file).rsplit('.', 1)[0]}_{streamflow_column}.geojson")

                    # Always regenerate the FIST GeoJSON so reruns pick up CRS and
                    # geometry fixes instead of silently reusing a stale output file.
                    if os.path.exists(GeoJSON_File):
                        os.remove(GeoJSON_File)
                    
                    tasks.append((
                        domain.vdt,
                        domain.stream_raster,
                        GeoJSON_File,
                        OutProjection,
                        domain.stream_geometry,
                        stream_id_field,
                        downstream_id_field,
                        SEED_File,
                        True, 
                        None,
                        streamflow_forecast_filtered_df
                    ))

    def get_consequences_tasks(self):
        tasks = []
        for domain in self.domains:
            for depth_file in domain.depth_maps:
                # Convert the depth raster to WGS84 coordinate system before processing
                depth_file_wgs84 = depth_file.replace('.tif', '_WGS84.tif')
                LOG.info(f"Converting depth raster to WGS84: {depth_file} -> {depth_file_wgs84}")
                try:
                    warp_options = gdal.WarpOptions(
                        format="GTiff",
                        dstSRS="EPSG:4326",
                        resampleAlg=gdal.GRA_Bilinear,
                        creationOptions=["COMPRESS=DEFLATE"],
                    )
                    ds = gdal.Warp(depth_file_wgs84, depth_file, options=warp_options)
                    if ds is None:
                        LOG.error(f"Failed to convert {depth_file} to WGS84, using original file")
                        depth_file_wgs84 = depth_file
                    else:
                        ds = None
                        LOG.info(f"Successfully converted to WGS84: {depth_file_wgs84}")
                except Exception as e:
                    LOG.error(f"Error converting to WGS84: {e}, using original file")
                    depth_file_wgs84 = depth_file
                
                Forecast_Flood_Depth_Raster_Name = os.path.basename(depth_file_wgs84)
                Consequences_JSON_File = Forecast_Flood_Depth_Raster_Name.replace('.tif','_consequences.json') 
                Consequences_JSON_Path = os.path.join(domain.folder_structure.Consequences_Folder, Consequences_JSON_File)
                Consequences_Output_GPKG_File = Consequences_JSON_File.replace('.json','.gpkg')
                LOG.info(f"Creating consequences file {Consequences_JSON_Path}")
                Create_Go_Consequence_GeoJSON(Consequences_JSON_Path, Forecast_Flood_Depth_Raster_Name, Consequences_Output_GPKG_File)
                # run the go-consequences Docker container
                source_dir = os.path.join(domain.folder_structure.top_level_folder, domain.name)
                source_dir = Path(source_dir).resolve().as_posix()

                docker_command = [
                    "docker",
                    "run",
                    "--rm",
                    "--mount",
                    f"type=bind,source={source_dir},target=/data",
                    "go-consequences:latest",
                    f"/data/Consequences/{Consequences_JSON_File}",
                ]
                tasks.append(docker_command)

    
    def make_bathymetry(self,
                        disable_bathymetry: bool,
                        mapper: Mapper = Mapper.CURVE2FLOOD_KERNEL_WEIGHTED,
                        overwrite: bool = False,
                        quiet: bool = False,):
        bathy_tasks = self.get_bathy_arc_tasks(bathymetry_disabled=disable_bathymetry, overwrite=overwrite)
        if bathy_tasks:
            with self.timer('arc_bathy'):
                for _ in tqdm.tqdm(self.iter(partial(_run_arc, quiet=quiet), bathy_tasks), desc="Running ARC for Bathymetry", disable=quiet):
                    pass

        bathy_mapper_tasks = self.get_bathy_mapper_tasks(bathymetry_disabled=disable_bathymetry, overwrite=overwrite)
        if bathy_mapper_tasks:
            with self.timer('flood_bathy'):
                for _ in tqdm.tqdm(self.iter(partial(_run_mapper, mapper=mapper, quiet=quiet), bathy_mapper_tasks), desc="Running Mapper for Bathymetry", disable=quiet):
                    pass

    def make_floodmaps(self, 
                       mapper: Mapper = Mapper.CURVE2FLOOD_KERNEL_WEIGHTED,
                       overwrite=False, 
                       quiet=False):
        floodmap_tasks = self.get_floodmap_tasks(overwrite=overwrite)
        if floodmap_tasks:
            with self.timer('flood'):
                for _ in tqdm.tqdm(self.iter(partial(_run_mapper, mapper=mapper, quiet=quiet), floodmap_tasks), desc="Running Mapper for Floodmaps", disable=quiet):
                    pass

    def make_FIST_inputs(self,
                          floodmap_mode: FloodMapMode,
                          streamflow_source: StreamflowSource,
                          stream_id_field: str,
                          downstream_id_field: str,
                          forensic_forecast_date: str | None,
                          forensic_forecast_hour: str | None,
                          forecastdate: str | None,
                          forecasthour: str | None,
                          overwrite=False,):
          FIST_tasks = self.get_FIST_tasks(
                floodmap_mode=floodmap_mode,
                streamflow_source=streamflow_source,
                stream_id_field=stream_id_field,
                downstream_id_field=downstream_id_field,
                forensic_forecast_date=forensic_forecast_date,
                forensic_forecast_hour=forensic_forecast_hour,
                forecastdate=forecastdate,
                forecasthour=forecasthour,
                overwrite=overwrite
          )
          if FIST_tasks:
                with self.timer('geojson_fist'):
                 for _ in tqdm.tqdm(self.iter(_run_fist, FIST_tasks), desc="Running Mapper for FIST inputs", disable=quiet):
                      pass
                 
    def run_go_consequences(self,
                            quiet=False):
        LOG.info("Creating the Go-Consequences JSON file and running the Go-Consequences Docker container...")

        consequences_tasks = self.get_consequences_tasks()
        if consequences_tasks:
            with self.timer('go_consequences'):
                for _ in tqdm.tqdm(self.iter(partial(subprocess.run, check=True), consequences_tasks), desc="Running Go-Consequences Docker container", disable=quiet):
                    pass
