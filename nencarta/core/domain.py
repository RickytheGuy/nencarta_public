import os
import io
import sys
import json
import time
import urllib
import warnings
import subprocess
from contextlib import nullcontext
from functools import cache
from typing import Literal
try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

os.environ["AWS_NO_SIGN_REQUEST"] = "YES"
os.environ["AWS_S3_ENDPOINT"] = "s3.amazonaws.com"
os.environ["KMP_WARNINGS"] = "0"

import yaml
import requests
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
from numba import njit
from pyproj import CRS
from osgeo import gdal, ogr
from shapely.geometry import box
from propcache import cached_property

from arc import Arc
from curve2flood import Curve2Flood_MainFunction

from .raster import Raster
from .vector import Vector
from nencarta.logger import LOG
from nencarta import DEM_Cleaner
from nencarta import Hydroterrain_Processing
from nencarta.core.configs import NencartaConfig
from nencarta.workspace import FolderStructure
from nencarta.core.enumerations import Mapper, StreamflowSource
from nencarta import Download_Process_ForecastData as ForecastFlows
from nencarta._constants import (
    GEOGLOWS_DAILY_URL, GEOGLOWS_FDC_URL, GEOGLOWS_RETURN_PERIODS_URL, NWM_RP_URL, CACHED_ESA_GRID
)

@cache
def get_rp_ds():
    """ This is faster for multiprocessing contexts since the dataset is only loaded once per process."""
    return xr.open_zarr(GEOGLOWS_RETURN_PERIODS_URL, storage_options={'anon': True})

@cache
def get_fdc_ds():
    return xr.open_zarr(GEOGLOWS_FDC_URL, storage_options={'anon': True})


@cache
def get_daily_ds():
    return xr.open_zarr(GEOGLOWS_DAILY_URL, storage_options={'anon': True})


def ignore_if_dead(method):
    def wrapper(self, *args, **kwargs):
        if self.dead:
            return self
        return method(self, *args, **kwargs)
    return wrapper

def _apply_buffer(bbox, distance):
    minx, miny, maxx, maxy = bbox
    return (
        minx - distance,
        miny - distance,
        maxx + distance,
        maxy + distance,
    )

def _build_options(bbox, vrt, xres=None, yres=None):
    option_dict = dict(
        resampleAlg="nearest",
        outputBounds=bbox,
    )

    if xres is not None:
        option_dict["xRes"] = xres
        option_dict["yRes"] = yres
        option_dict["targetAlignedPixels"] = True

    if vrt:
        option_dict["outputSRS"] = "EPSG:4326"
        return gdal.BuildVRT, gdal.BuildVRTOptions(**option_dict)

    option_dict["dstSRS"] = "EPSG:4326"
    return gdal.Warp, gdal.WarpOptions(**option_dict)

@njit(cache=True)
def _clean_stream_raster_array(array: np.ndarray) -> tuple[np.ndarray, int, int]:
    mask = array > 0
    array = np.where(mask, array, 0)
    (RR,CC) = np.where(mask)
    num_nonzero = len(RR)
    first_pass = 0
    second_pass = 0
    
    for filterpass in range(2):
        #First pass is just to get rid of single cells hanging out not doing anything
        p_count = 0
        p_percent = (num_nonzero+1)/100.0
        for x in range(num_nonzero):
            if x>=p_count*p_percent:
                p_count = p_count + 1
            r=RR[x]
            c=CC[x]
            V = array[r,c]
            if V>0:
                #Left and Right cells are zeros
                if array[r,c+1]==0 and array[r,c-1]==0:
                    #The bottom cells are all zeros as well, but there is a cell directly above that is legit
                    if (array[r+1,c-1]+array[r+1,c]+array[r+1,c+1])==0 and array[r-1,c]>0:
                        array[r,c] = 0
                        first_pass += 1
                    #The top cells are all zeros as well, but there is a cell directly below that is legit
                    elif (array[r-1,c-1]+array[r-1,c]+array[r-1,c+1])==0 and array[r+1,c]>0:
                        array[r,c] = 0
                        first_pass += 1
                #top and bottom cells are zeros
                if array[r,c]>0 and array[r+1,c]==0 and array[r-1,c]==0:
                    #All cells on the right are zero, but there is a cell to the left that is legit
                    if (array[r+1,c+1]+array[r,c+1]+array[r-1,c+1])==0 and array[r,c-1]>0:
                        array[r,c] = 0
                        first_pass += 1
                    elif (array[r+1,c-1]+array[r,c-1]+array[r-1,c-1])==0 and array[r,c+1]>0:
                        array[r,c] = 0
                        first_pass += 1
        
        #This pass is to remove all the redundant cells
        p_count = 0
        p_percent = (num_nonzero+1)/100.0
        for x in range(num_nonzero):
            if x>=p_count*p_percent:
                p_count = p_count + 1
            r=RR[x]
            c=CC[x]
            V = array[r,c]
            if V>0:
                if array[r+1,c]==V and (array[r+1,c+1]==V or array[r+1,c-1]==V):
                    if sum(array[r+1,c-1:c+2])==0:
                        array[r+1,c] = 0
                        second_pass += 1
                elif array[r-1,c]==V and (array[r-1,c+1]==V or array[r-1,c-1]==V):
                    if sum(array[r-1,c-1:c+2])==0:
                        array[r-1,c] = 0
                        second_pass += 1
                elif array[r,c+1]==V and (array[r+1,c+1]==V or array[r-1,c+1]==V):
                    if sum(array[r-1:r+1,c+2])==0:
                        array[r,c+1] = 0
                        second_pass += 1
                elif array[r,c-1]==V and (array[r+1,c-1]==V or array[r-1,c-1]==V):
                    if sum(array[r-1:r+1,c-2])==0:
                            array[r,c-1] = 0
                            second_pass += 1

    return array, first_pass, second_pass

def clean_stream_raster(stream_ds: gdal.Dataset, num_passes: int = 2) -> None:
    """
    This function comes from Mike Follum's ARC at https://github.com/MikeFHS/automated-rating-curve
    """
    assert num_passes > 0, "num_passes must be greater than 0"
    
    # Get stream raster
    shape = (stream_ds.RasterYSize + 2, stream_ds.RasterXSize + 2)
    array = np.zeros(shape, dtype=np.int64)
    array[1:-1, 1:-1] = stream_ds.ReadAsArray()

    array, first_pass, second_pass = _clean_stream_raster_array(array)
    LOG.info(f'First Pass - Removed {first_pass} cells.')
    LOG.info(f'Second Pass - Removed {second_pass} cells.')
    
    # Write the cleaned array to the raster
    stream_ds.WriteArray(array[1:-1, 1:-1])

def get_geoglows_rp(river_ids: list[int]) -> pd.DataFrame:
    rp_ds = get_rp_ds().sel(river_id=river_ids)[['logpearson3', 'gumbel']]
        
    # Convert Xarray to Dask DataFrame and pivot
    rp_df = rp_ds.to_dataframe().reset_index()

    # find the maximum between the gumbel and logpearson3 return periods and label this new column 'return_period_flow'
    rp_df['return_period_flow'] = rp_df[['gumbel', 'logpearson3']].mean(axis=1)

    # keep just the column 'return_period_flow'
    rp_df = rp_df[['river_id', 'return_period', 'return_period_flow']]

    # Convert 'return_period' to category dtype
    rp_df['return_period'] = rp_df['return_period'].astype('category')
    
    # Pivot the table
    rp_df = rp_df.pivot_table(index='river_id', columns='return_period', values='return_period_flow', aggfunc='mean', observed=False)

    # Rename columns to indicate return periods
    rp_df = rp_df.rename(columns={col: f'rp{int(col)}' for col in rp_df.columns})

    p_exceedances = np.arange(0, 105, 5, dtype=float)
    try:
        fdc_ds = get_fdc_ds().sel(p_exceed=p_exceedances, river_id=river_ids)

        # Convert Xarray to Dask DataFrame
        fdc_df = fdc_ds.to_dataframe().reset_index()

        fdc_df = fdc_df.pivot_table(
            index='river_id',
            columns='p_exceed',
            values='hourly_annual',
            aggfunc='mean'
        )
        fdc_df = fdc_df.rename(columns={p: f"p_exceed_{p}" for p in fdc_df.columns})
    except:
        LOG.warning("FDC data not available; falling back to daily data for FDC calculation.")
        # Load daily data from S3 using Dask
        # Convert to a list of integers
        dailyflow_ds = get_daily_ds().sel(river_id=river_ids)
        # Convert Xarray to Dask DataFrame
        daily_df = dailyflow_ds.to_dataframe().reset_index()

        # creating exceedance percentiles with the daily data
        quantiles = [1.0 - (p / 100.0) for p in p_exceedances]
        fdc_df = daily_df.groupby('river_id')['Q'].quantile(quantiles).unstack()
        fdc_df = fdc_df.rename(
            columns={q: f"p_exceed_{p}" for q, p in zip(quantiles, p_exceedances)}
        )

        # uniqify the index
        fdc_df = fdc_df[~fdc_df.index.duplicated(keep='first')]

    final_df = pd.concat([fdc_df, rp_df], axis=1)
    final_df['COMID'] = final_df.index

    # Reorder the DataFrame
    columns = ['COMID'] + [col for col in final_df.columns if col != 'COMID']
    final_df = final_df[columns]

    for col in ['p_exceed_0', 'rp100']:
        # I think this is a better way of buffering the maximum flow
        # Multiping by 1.5 seems to be a reasonable esimate of the maximum high flow, while adding 50 helps small rivers with tiny
        # return period 100 flows (close to 0)
        # Going too big means the VDT has bigger gaps to fill, which can lead to worse performance and less accurate rating curves
        final_df[f'{col}_premium'] = final_df[col] * 1.5 + 50

    final_df = final_df.round(3)
    return final_df

def get_nwm_rp(comids: list[int], nwm_api_key: str):
    if not nwm_api_key:
        raise ValueError("nwm_api_key is required for NWM return period requests.")

    header = {'x-api-key': nwm_api_key}
    params = {'comids': ','.join(map(str, comids)),
              'output_format': 'csv',
              'order_by_comid': False,}

    response = requests.get(NWM_RP_URL, params=params, headers=header, timeout=60)

    if response.status_code == 200:
        return_period_df = pd.read_csv(io.StringIO(response.text))
    else:
        raise requests.exceptions.HTTPError(response.text)
    
    return_period_df = return_period_df.set_index("feature_id")
    return_period_df.index.name = "river_id"
    return_period_df.columns = ['rp2', 'rp5', 'rp10', 'rp25', 'rp50', 'rp100']

    # Add derived flows directly to rp_df without dropping anything
    return_period_df["rp100_premium"] = (return_period_df["rp100"] * 1.5) + 50

    # Reorder columns so the return period fields come first
    cols = [col for col in return_period_df.columns if col.startswith("rp")]
    return_period_df = return_period_df[cols]

    return_period_df['COMID'] = return_period_df.index

    # Reorder the DataFrame
    columns = ['COMID'] + [col for col in return_period_df.columns if col != 'COMID']
    return_period_df = return_period_df[columns]

    return return_period_df

def filter_streams_with_lake_json(lake_filter_json: str, stream_df: gpd.GeoDataFrame, rivid_field: str) -> gpd.GeoDataFrame:
    # First let's remove the stream reaches that are in the stream_ids_in_lake_list
    # filter out the streams that are in the stream_ids_in_lake_list by using the "LINKNO values in stream_ids_in_lake_list"
    with open(lake_filter_json, 'r') as f:
        lake_filter: dict[str, dict[str, list[int]]] = json.load(f)
        stream_ids_in_lake_list = []
        for _k, v in lake_filter.items():
            inside = v.get("inside", [])
            for x in inside:
                if x is not None:
                    stream_ids_in_lake_list.append(x)
    return stream_df[~stream_df[rivid_field].isin(stream_ids_in_lake_list)]

def filter_streams_by_stream_order(stream_df: gpd.GeoDataFrame, strm_order_field: str, strm_order_low: float | None, strm_order_high: float | None) -> gpd.GeoDataFrame:
    if strm_order_field not in stream_df.columns:
        LOG.warning(f"StrmOrder_Field '{strm_order_field}' not found in stream shapefile; skipping stream order filter.")
        return stream_df

    stream_df[strm_order_field] = pd.to_numeric(stream_df[strm_order_field], errors="coerce")
    if strm_order_low is None:
        strm_order_low = stream_df[strm_order_field].min()
    if strm_order_high is None:
        strm_order_high = stream_df[strm_order_field].max()

    return stream_df[stream_df[strm_order_field].between(strm_order_low, strm_order_high)]

def download_grid_with_retry(url, max_retries=10, wait_seconds=10):
    attempt = 0
    while attempt < max_retries:
        try:
            grid = gpd.read_file(url).to_crs(epsg=4326)
            LOG.info(f"Successfully downloaded grid on attempt {attempt + 1}")
            return grid
        except urllib.error.URLError as e:
            attempt += 1
            LOG.warning(f"Attempt {attempt} failed: {e}")
            time.sleep(wait_seconds)
        except Exception as e:
            attempt += 1
            LOG.warning(f"Unexpected error on attempt {attempt}: {e}")
            time.sleep(wait_seconds)

    raise RuntimeError(f"Failed to download grid after {max_retries} attempts.")

def get_esa_tiles(bbox):
    if os.path.isfile(CACHED_ESA_GRID):
        try:
            return set(gpd.read_parquet(CACHED_ESA_GRID, bbox=bbox)['ll_tile'])
        except ValueError as e:
            if "bbox" in str(e):
                grid = gpd.read_parquet(CACHED_ESA_GRID)
                grid.to_parquet(CACHED_ESA_GRID, write_covering_bbox=True)
                return set(grid[grid.intersects(box(*bbox))]['ll_tile'])

            raise e
    
    url = "https://esa-worldcover.s3.eu-central-1.amazonaws.com/esa_worldcover_grid.geojson"
    grid = download_grid_with_retry(url)
    grid.to_parquet(CACHED_ESA_GRID, write_covering_bbox=True)
    return set(grid[grid.intersects(box(*bbox))]['ll_tile'])

class Domain():
    def __init__(self, folder_structure: FolderStructure):
        self.folder_structure = folder_structure
        self.dem = None
        self._dem_basename = None
        self.source_dems: list[str] = []
        self.stream_geometry = None
        self.source_stream_geometry: list[str] = []
        self.stream_attribute = None
        self.land_cover = None
        self.stream_raster = None
        self.reanalysis_flow_file = None
        self.arc_config = None
        self.baseflow_file = None
        self.baseflow_floodmap = None
        self.flood_flow_files: list[str] = []

        # Floodmap outputs
        self.floodmap_configs: list[str] = []
        self.flood_maps: list[str] = []
        self.floodmap_vectors: list[str] = []
        self.velocity_maps: list[str] = []
        self.wse_maps: list[str] = []
        self.depth_maps: list[str] = []

        # ARC outputs
        self.vdt = None
        self.curve_file = None
        self.ap_file = None
        self.arc_bathy = None

        # FLDPLN inputs
        self.filled_dem = None
        self.flow_direction = None
        self.flow_accumulation = None
        self.new_stream_geometry = None
        self.new_catchment = None
        self.new_stream_geometry_matched = None

        # Clean DEM related
        self.clean_dem_config = None
        self.vdt_for_clean_dem = None
        self.curve_file_for_clean_dem = None
        self.initial_floodmap_for_clean_dem = None
        self.floodmap_gpkg_for_clean_dem = None
        self.cleaned_dem = None
        self.water_mask = None

        # Bathymetry related
        self.arc_bathymetry = None
        self.dem_with_bathymetry = None

        self.dead = False

    @property
    def name(self) -> str:
        return os.path.basename(self.folder_structure.top_level_folder)

    def assign_source_dem(self, dem: str | list[str]) -> Self:
        if isinstance(dem, str):
            self.source_dems.append(dem)
        elif isinstance(dem, list):
            self.source_dems.extend(dem)
        else:
            raise ValueError("DEM must be a string or a list of strings.")
        
        return self

    def __getstate__(self):
        state = self.__dict__.copy()
        # Don't pickle the raster objects
        if 'dem_raster' in state:
            del state['dem_raster']
        return state

    def get_priority(self):
        if self.dead:
            return 0
        
        if not self.dem:
            raise ValueError("DEM must be assigned before calculating priority.")
        
        return -os.path.getsize(self.dem)

    @cached_property
    def dem_raster(self) -> Raster:
        if not self.dem:
            raise ValueError("DEM must be assigned before accessing the raster.")
        return Raster(self.dem)
    
    @property
    @ignore_if_dead
    def dem_basename(self) -> str:
        if not self._dem_basename:
            raise ValueError("DEM must be assigned before accessing the basename.")
        
        return self._dem_basename

    def assign_dem(self, 
                   dem: str = None, 
                   bbox: tuple = None, 
                   buffer: bool = False,
                   vrt: bool = False, 
                   buffer_distance: float = 0.05, 
                   raise_error_if_no_dems: bool = True,
                   overwrite: bool = False) -> Self:
        if buffer and not self.source_dems:
            raise ValueError(
                "Buffering requested but no source DEMs assigned."
            )

        if dem:
            basename = os.path.splitext(os.path.basename(dem))[0]
        else:
            basename = f"dem_{bbox[0]}_{bbox[1]}_{bbox[2]}_{bbox[3]}"
        output_dem = os.path.join(
            self.folder_structure.dem_dir,
            f"{'buffered_' if buffer else ''}{basename}.{'vrt' if vrt else 'tif'}"
        )

        self.dem = output_dem
        self._dem_basename = basename

        if os.path.exists(output_dem) and not overwrite:
            return self

        os.makedirs(os.path.dirname(output_dem), exist_ok=True)

        if dem and bbox:
            target_bbox = bbox

        elif dem and buffer:
            raster = Raster(dem)
            target_bbox = raster.epsg_4326_bbox

        elif bbox:
            target_bbox = bbox

        elif dem:
            if vrt:
                gdal.BuildVRT(output_dem, dem)
            else:
                self._dem_basename = os.path.splitext(os.path.basename(dem))[0]
                self.dem = dem

            return self

        else:
            raise ValueError(
                "Must specify either a DEM or a bounding box."
            )

        if buffer:
            target_bbox = _apply_buffer(target_bbox, buffer_distance)

        if dem:
            candidates = self.source_dems + [dem] if self.source_dems else [dem]

            surrounding_dems = Raster.get_rasters_in_extent(
                target_bbox,
                candidates,
            )

            assert dem in surrounding_dems

        else:
            surrounding_dems = Raster.get_rasters_in_extent(
                target_bbox,
                self.source_dems,
            )

            if not surrounding_dems:
                if raise_error_if_no_dems:
                    raise ValueError(
                        "No source DEMs intersect the specified bounding box."
                    )
                
                self.dead = True
                return self

        xres = yres = None

        if dem:
            raster = Raster(dem)
            xres, yres = raster.resolution
        else:
            # Use the resolution of the first intersecting DEM as the target resolution
            # This avoids little gaps when stitching rasters files together
            raster = Raster(surrounding_dems[0])
            xres, yres = raster.resolution

        builder, options = _build_options(
            target_bbox,
            vrt=vrt,
            xres=xres,
            yres=yres,
        )

        builder(output_dem, surrounding_dems, options=options)

        return self
    
    def assign_flood_flow_files(self, flow_files: list[str]) -> Self:
        self.flood_flow_files.extend(flow_files)
        return self
    
    @ignore_if_dead
    def validate_dem_coordinate_system(self):
        if not self.dem:
            raise ValueError("DEM must be assigned before validating coordinate system.")
        
        #Check the coordinate system of the rasters and if they are not in meters or degrees, end with log an error message and stop processing
        unit_aliases = {
            'meter': 'meter',
            'meters': 'meter',
            'metre': 'meter',
            'metres': 'meter',
            'degree': 'degree',
            'degrees': 'degree'
        }
        raster_projections = {
            'DEM': self.dem_raster.projection,
        }
        for raster_name, raster_projection in raster_projections.items():
            try:
                raster_crs = CRS.from_wkt(raster_projection)
            except Exception as ex:
                LOG.error(f'Unable to parse CRS for {raster_name} raster: {ex}')
                return False

            axis_units = {(axis.unit_name or '').strip().lower() for axis in raster_crs.axis_info if axis is not None}
            axis_units.discard('')
            if not axis_units:
                LOG.error(f'Unable to determine CRS units for {raster_name} raster.')
                return False

            invalid_units = [u for u in sorted(axis_units) if unit_aliases.get(u) not in {'meter', 'degree'}]
            if invalid_units:
                LOG.error(f'{raster_name} raster CRS units are not meters or degrees: {", ".join(invalid_units)}')
                return False
            
        return True
    
    @ignore_if_dead
    def make_stream_geometry(self,
                                 name: str,
                                 stream_geometry: str | list[str],
                                 raise_error_if_no_streams: bool = True,
                                 overwrite: bool = False) -> Self:
        if not self.dem:
            raise ValueError("DEM must be assigned before generating stream raster.")
        
        self.stream_geometry = os.path.join(self.folder_structure.strm_folder, name)
        os.makedirs(os.path.dirname(self.stream_geometry), exist_ok=True)

        if not os.path.exists(self.stream_geometry) or overwrite:
            bbox = self.dem_raster.epsg_4326_bbox

            streamlines = Vector.get_streamlines_in_extent(bbox, stream_geometry if isinstance(stream_geometry, list) else [stream_geometry])
            if not streamlines:
                if raise_error_if_no_streams:
                    raise ValueError("No RFS stream geometry files intersect the DEM extent.")
                
                self.dead = True
                return self

        if not os.path.exists(self.stream_geometry) or overwrite:
            if len(streamlines) == 1:
                gdf = Vector(streamlines[0]).to_geopandas(bbox_epsg_4326=bbox)
            else:
                gdf = pd.concat([Vector(path).to_geopandas(bbox_epsg_4326=bbox) for path in streamlines], ignore_index=True, copy=False)

            gdf = gdf[~gdf.geometry.isna()]
            # Legacy NenCarta expects both LINKNO and COMID to exist for downstream
            # tools/configs, even when the source stream dataset only contains one.
            if "LINKNO" not in gdf.columns and "COMID" in gdf.columns:
                gdf["LINKNO"] = gdf["COMID"]
            if "COMID" not in gdf.columns and "LINKNO" in gdf.columns:
                gdf["COMID"] = gdf["LINKNO"]

            if gdf.empty:
                if raise_error_if_no_streams:
                    raise ValueError("No stream geometries intersect the DEM extent.")
                
                self.dead = True
                return self
            
            kwargs = {}
            if self.stream_geometry.endswith(('.pq', '.parquet')):
                kwargs['compression'] = 'brotli'
                kwargs['write_covering_bbox'] = True

            Vector.save_any_geom(gdf, self.stream_geometry, **kwargs)

        return self

    @ignore_if_dead
    def make_fldpln_inputs(self, 
                               new_strm_threshold_km2: float, 
                               stream_id_field: str,
                               ds_stream_id_field: str,
                               strm_order_field: str,
                               min_match_score: float,
                               overwrite: bool = False) -> Self:
        if not self.dem or not self.stream_geometry:
            raise ValueError("DEM and stream geometry must be assigned before generating FLDPLN inputs.")
        
        self.filled_dem = os.path.join(self.folder_structure.Flow_Direction_Folder, f"{self.dem_basename}_filled.tif")
        self.flow_direction = os.path.join(self.folder_structure.Flow_Direction_Folder, f"{self.dem_basename}_flowdir.tif")
        self.flow_accumulation = os.path.join(self.folder_structure.Flow_Direction_Folder, f"{self.dem_basename}_flowacc.tif")
        self.new_stream_geometry = os.path.join(self.folder_structure.Flow_Direction_Folder, f"{self.dem_basename}_flowlines.gpkg")
        self.new_catchment = os.path.join(self.folder_structure.Flow_Direction_Folder, f"{self.dem_basename}_catchments.gpkg")
        self.new_stream_geometry_matched = os.path.join(self.folder_structure.strm_folder, f"{self.dem_basename}_flowlines_matched.gpkg")

        
        if not overwrite and os.path.exists(self.new_stream_geometry_matched) and os.path.exists(self.filled_dem):
            return self
        
        os.makedirs(self.folder_structure.Flow_Direction_Folder, exist_ok=True)
        os.makedirs(self.folder_structure.strm_folder, exist_ok=True)

        if overwrite:
            for vector_path in (
                self.new_stream_geometry,
                self.new_catchment,
                self.new_stream_geometry_matched,
            ):
                Vector.delete_any_geom(vector_path)

        Hydroterrain_Processing.create_flow_direction_and_flow_accumulation_raster(
            self.dem, self.filled_dem, self.flow_direction, self.flow_accumulation
        )
        Hydroterrain_Processing.create_catchments_and_flowlines_with_flow_direction_and_accumulation(
            self.flow_direction,
            self.flow_accumulation,
            self.folder_structure.Flow_Direction_Folder,
            new_strm_threshold_km2,
            self.new_stream_geometry,
            self.new_catchment
        )
        Hydroterrain_Processing.match_new_streams_to_old_streams(
            self.new_stream_geometry,
            self.stream_geometry,
            self.new_stream_geometry_matched,
            stream_id_field,
            ds_stream_id_field,
            strm_order_field,
            min_match_score = min_match_score,
            require_overlap = False,
            remove_detached_upstream = True,
            connectivity_tolerance_m = 30.0,
            buffer_distance_m =1000.0
        )

        self.dem = self.filled_dem
        del self.dem_raster # Clear cached property so it will be reloaded from the new filled DEM path on next access
        self.stream_geometry = self.new_stream_geometry_matched
    
    @ignore_if_dead
    def make_stream_raster(self, 
                               name: str,
                               river_id: str = 'LINKNO', 
                               overwrite: bool = False) -> Self:
        if not self.stream_geometry:
            raise ValueError("Stream geometry must be assigned before generating stream raster.")
        
        self.stream_raster = os.path.join(self.folder_structure.strm_folder, name)

        if os.path.exists(self.stream_raster) and not overwrite:
            return self
        
        os.makedirs(self.folder_structure.strm_folder, exist_ok=True)

        warnings.filterwarnings(
            "once",
            message="Failed to fetch spatial reference on layer.*",
            category=RuntimeWarning,
        )

        dem_raster: Raster = self.dem_raster
        driver: gdal.Driver = gdal.GetDriverByName('GTiff')

        if os.path.exists(self.stream_raster):
            try:
                delete_failed = driver.Delete(self.stream_raster) != 0
            except (PermissionError, RuntimeError, OSError):
                delete_failed = True

            if delete_failed:
                raise PermissionError(
                    f"Unable to overwrite stream raster because it is locked: {self.stream_raster}"
                )

        source_ds = gdal.OpenEx(self.stream_geometry, gdal.OF_VECTOR)
        if source_ds is None:
            raise RuntimeError(f"Unable to open stream geometry: {self.stream_geometry}")

        layer: ogr.Layer = source_ds.GetLayer(0)
        layer_name = layer.GetName() if layer is not None else None

        stream_ds: gdal.Dataset = driver.Create(self.stream_raster, dem_raster.ncols, dem_raster.nrows, 1, gdal.GDT_Int32, ['COMPRESS=DEFLATE', 'PREDICTOR=2'])
        stream_ds.SetGeoTransform(dem_raster.geotransform)
        stream_ds.SetProjection(dem_raster.projection)
        options = gdal.RasterizeOptions(
            attribute=river_id,
            layers=[layer_name] if layer_name and self.stream_geometry.lower().endswith(".gpkg") else None,
        )
        gdal.Rasterize(stream_ds, source_ds, options=options)
        source_ds = None

        # Clean the raster
        clean_stream_raster(stream_ds)
        stream_ds = None

        return self
    
    @ignore_if_dead
    def make_land_cover(self, 
                        name: str,
                        land_cover_cache: list[str] = None, 
                        year: Literal[2020, 2021] = 2021, 
                        overwrite: bool = False) -> Self:
        if not self.dem:
            raise ValueError("DEM must be assigned before generating land cover.")
        
        self.land_cover = os.path.join(self.folder_structure.land_folder, name)
        if os.path.exists(self.land_cover) and not overwrite:
            LOG.info(f"Land cover file already exists at {self.land_cover}, skipping generation.")
            return self
        
        os.makedirs(os.path.dirname(self.land_cover), exist_ok=True)

        bounds = self.dem_raster.epsg_4326_bbox
        tiles = get_esa_tiles(bbox=bounds)

        version_by_year = {2020: "v100", 2021: "v200"}  # extend as needed for newer years
        if year not in version_by_year:
            raise ValueError(f"Year {year} not supported. Available: {list(version_by_year)}")
        version = version_by_year[year]

        landcover_files: list[str] = []
        temp_downloaded: list[str] = []
        cached_files = {}
        if land_cover_cache:
            cached_files = {os.path.splitext(os.path.basename(path))[0]: path for path in land_cover_cache}
        for tile in tiles:
            if tile in cached_files:
                landcover_files.append(cached_files[tile])
            else:
                landcover_files.append(f"/vsis3/esa-worldcover/{version}/{year}/map/ESA_WorldCover_10m_{year}_{version}_{tile}_Map.tif")

        LOG.info(f"Creating {self.land_cover}")
        if not landcover_files:
            out_ds: gdal.Dataset = gdal.GetDriverByName('GTiff').Create(self.land_cover, self.dem_raster.ncols, self.dem_raster.nrows, 1, gdal.GDT_Byte, {'COMPRESS': 'DEFLATE', 'PREDICTOR': '2'})
            out_ds.SetGeoTransform(self.dem_raster.geotransform)
            out_ds.SetProjection(self.dem_raster.projection)
            LOG.warning("No land cover files found for the specified DEM extent. Generating a fake land cover file filled with 10 (trees).")
            # Let's make a fake landcover file for arc to use. 
            # Fill it with 10, since the areas that don't have it tend to be tropical islands (10 is trees)
            out_ds.GetRasterBand(1).Fill(10)
            return
        
        proj = self.dem_raster.projection
        xres, yres = self.dem_raster.resolution
        if name.endswith('.vrt'):
            options = gdal.BuildVRTOptions(outputBounds=bounds,
                                outputSRS=proj,
                                xRes=xres,
                                yRes=yres,
                                resampleAlg='mode')
            gdal.BuildVRT(self.land_cover, landcover_files, options=options)
        else:
            options = gdal.WarpOptions(format='GTiff',
                                outputBounds=bounds,
                                outputBoundsSRS=proj,
                                dstSRS=proj,
                                xRes=xres,
                                yRes=yres,
                                resampleAlg='mode',
                                creationOptions=['COMPRESS=DEFLATE', 'PREDICTOR=2'])
            gdal.Warp(self.land_cover, landcover_files, options=options)

        return self

    @ignore_if_dead
    def make_base_max_flows(self, 
                                rivid_field: str,
                                filename: str, 
                                flow_source: StreamflowSource,
                                highflow_field: str,
                                strm_order_field: str = None,
                                strm_order_low: int = None,
                                strm_order_high: int = None,
                                baseflow_field: str = None,
                                baseflow_threshold: float = None,
                                lake_filter_json: str = None,
                                nwm_api_key: str = None,
                                overwrite: bool = False) -> Self:
        """
        This function generates a CSV file containing base and maximum flow values for each stream segment in the domain, based on the stream geometry and precomputed flow datasets. The flow values are derived from both the Flow Duration Curve (FDC) and Return Period (RP) datasets, which are accessed via Dask arrays for efficient computation. The resulting CSV file includes columns for various return periods and exceedance probabilities, as well as "premium" flow values calculated as 1.5 times the base flow plus 50.
        This is inspired by nencarta's equivalent function.
        """
        if not self.stream_geometry:
            raise ValueError("Stream geometry must be generated before generating base/max flows.")
        
        self.reanalysis_flow_file = os.path.join(self.folder_structure.FLOW_Folder, filename)
        if os.path.exists(self.reanalysis_flow_file) and not overwrite:
            return self

        os.makedirs(os.path.dirname(self.reanalysis_flow_file), exist_ok=True)
        stream_df = Vector(self.stream_geometry).to_geopandas()

        if lake_filter_json:
            stream_df = filter_streams_with_lake_json(lake_filter_json, stream_df, rivid_field)

        # if the StrmOrder_Field and StrmOrder_Lower or StrmOrder_Upper are not None use these to filter the StrmShp_gdf
        if strm_order_field and (strm_order_low is not None or strm_order_high is not None):
            stream_df = filter_streams_by_stream_order(stream_df, strm_order_field, strm_order_low, strm_order_high)

        river_ids = stream_df[rivid_field].astype(int).unique()

        if len(river_ids) == 0:
            LOG.warning("No stream segments remain after filtering; cannot generate base/max flow file.")
            self.dead = True
            return self

        if flow_source == StreamflowSource.GEOGLOWS:
            final_df = get_geoglows_rp(river_ids)
        elif flow_source.is_nwm():
            nwm_api_key = nwm_api_key or os.getenv("NWM_API_KEY")
            if not nwm_api_key:
                raise ValueError("NWM_API_KEY environment variable must be set for NWM flow retrieval.")
            final_df = get_nwm_rp(river_ids, nwm_api_key)
        else:
            raise ValueError("Invalid flow_source specified. Must be 'geoglows' or 'nwm'.")
    
        # break the code if the dataframe is empty or if the streamflow is all 0
        if final_df.empty or final_df[highflow_field].values.mean() <= 0:
            LOG.warning(f"Results for {self.dem_basename} are not possible because we don't have streamflow estimates...")
            self.dead = True
            return

        LOG.info(final_df)

        if baseflow_threshold is not None:
            if baseflow_field not in final_df.columns:
                LOG.warning(
                    f"baseflow_threshold was provided ({baseflow_threshold}), but baseflow field "
                    f"'{baseflow_field}' was not found in streamflow data. Skipping baseflow threshold filter."
                )
            else:
                final_df_before_filter_count = len(final_df)
                final_df = final_df[final_df[baseflow_field] >= baseflow_threshold]
                if final_df.empty:
                    LOG.warning(
                        f"All streams were removed by baseflow_threshold={baseflow_threshold} "
                        f"using field '{baseflow_field}'."
                    )
                    self.dead = True
                    return
                
                LOG.info(
                    f"Filtered out {final_df_before_filter_count - len(final_df)} streams below baseflow threshold of {baseflow_threshold} using field '{baseflow_field}'."
                )
                stream_df = stream_df[stream_df[rivid_field].isin(final_df['COMID'])]

        # Overwrite the stream geometry with the filtered one, so that downstream processes only have access to the streams for which we have flow data
        kwargs = {'index': False}
        if self.stream_geometry.endswith(('.pq', '.parquet')):
            kwargs['compression'] = 'brotli'
            kwargs['write_covering_bbox'] = True

        Vector.save_any_geom(stream_df, self.stream_geometry, **kwargs)

        if self.reanalysis_flow_file.endswith('.parquet'):
            final_df.round(3).to_parquet(self.reanalysis_flow_file, index=False)
        else:
            final_df.round(3).to_csv(self.reanalysis_flow_file, index=False)

        return self
    
    def _arc_inputs(self,
                    dem: str,
                    stream_raster: str,
                    land_cover: str,
                    mannings_n_file: str,
                    reanalysis_flow_file: str,
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

    def _fldpln_inputs(self):
        if not self.flow_direction:
            raise ValueError("Flow direction must be generated before defining FLDPLN inputs.")
        
        return {
            "# FLDPLN_Inputs": "",
            "Flow_Direction_File": self.flow_direction,
        }
    
    def write_config(self, config_path: str, params: dict):
        with open(config_path, 'w') as f:
            if config_path.endswith(('.yaml', '.yml')):
                    yaml.dump(params, f)
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

    @ignore_if_dead
    def define_configs_for_dem_cleaning(self,
                                        config_name: str,
                                        vdt_name: str,
                                        curve_name: str, 
                                        floodmap_name: str,
                                        floodmap_vector_name: str,
                                        cleaned_dem_name: str,
                                        mapper: Mapper = Mapper.CURVE2FLOOD_KERNEL_WEIGHTED,
                                        baseflow: str = 'p_exceed_50',
                                        maxflow: str = 'rp100_premium',
                                        reach_average_curve_file: bool = False,
                                        bathy_use_banks: bool = False,
                                        find_banks_from_land_cover: bool = False,
                                        flood_waterlc_and_strm_cells: bool = False,
                                        land_cover_water_value: int = 80,
                                        bathy_args: dict = {},
                                        use_specified_depths: bool = False,
                                        specified_depth: float = None,
                                        overwrite: bool = False
                                        ) -> Self:
        if not self.reanalysis_flow_file:
            raise ValueError("Base/max flow file must be generated before defining configs for DEM cleaning.")
        
        if not self.land_cover:
            raise ValueError("Land cover must be generated before defining configs for DEM cleaning.")

        if not self.baseflow_file:
            raise ValueError("Baseflow file must be generated before defining configs for DEM cleaning.")

        self.clean_dem_config = os.path.join(self.folder_structure.ARC_Folder, config_name)
        self.vdt_for_clean_dem = os.path.join(self.folder_structure.VDT_Folder, vdt_name)
        self.curve_file_for_clean_dem = os.path.join(self.folder_structure.VDT_Folder, curve_name)
        self.initial_floodmap_for_clean_dem = os.path.join(self.folder_structure.flood_folder, floodmap_name)
        self.floodmap_gpkg_for_clean_dem = os.path.join(self.folder_structure.flood_folder, floodmap_vector_name)
        self.cleaned_dem = os.path.join(self.folder_structure.dem_updated_folder, cleaned_dem_name)

        os.makedirs(self.folder_structure.ARC_Folder, exist_ok=True)
        os.makedirs(self.folder_structure.VDT_Folder, exist_ok=True)
        os.makedirs(self.folder_structure.flood_folder, exist_ok=True)
        os.makedirs(self.folder_structure.dem_updated_folder, exist_ok=True)

        if os.path.exists(self.clean_dem_config) and not overwrite:
            return self

        params = self._arc_inputs(
            dem=self.dem,
            stream_raster=self.stream_raster,
            land_cover=self.land_cover,
            mannings_n_file=self.folder_structure.mannings_text_file,
            reanalysis_flow_file=self.reanalysis_flow_file,
            river_id="COMID",
            highflow_field=maxflow,
            baseflow_field=baseflow,
            bathy_args=bathy_args,
            include_baseflow=True,
        )

        params["#VDT_Output_File_and_CurveFile"] = ""
        params["VDT_Database_NumIterations"] = 30
        params["Print_VDT_Database"] = self.vdt_for_clean_dem
        params["Print_Curve_File"] = self.curve_file_for_clean_dem
        params["Reach_Average_Curve_File"] = reach_average_curve_file

        params["#Mapper Input Data"] = ""
        params["StrmShp_File"] = self.stream_geometry
        params["Comid_Flow_File"] = self.baseflow_file
        params["mapper"] = mapper.value
        params["FS_ADJUST_FLOW_BY_FRACTION"] = 1.0
        params["Bathy_Use_Banks"] = bathy_use_banks
        if flood_waterlc_and_strm_cells:
            params["LAND_WaterValue"] = land_cover_water_value
        if find_banks_from_land_cover:
            params["FindBanksBasedOnLandCover"] = True
        if use_specified_depths:
            if mapper.is_curve2flood_fldpln_mapper():
                params.update(self._fldpln_inputs())
            params["OutFLD"] = self.initial_floodmap_for_clean_dem
            params["OutSHP"] = self.floodmap_gpkg_for_clean_dem
            params["FloodSpreader_SpecifyDepth"] = specified_depth

        self.write_config(self.clean_dem_config, params)
        return self

    @ignore_if_dead
    def define_arc_configs(self, 
                           config_name: str,
                           vdt_name: str,
                           floodmap_name: str,
                           floodmap_vector_name: str,
                           curve_file_name: str = None,
                           ap_file_name: str = None,
                           baseflow: str = 'p_exceed_50',
                           maxflow: str = 'rp100_premium',
                           reach_average_curve_file: bool = False,
                           mapper: Mapper = Mapper.CURVE2FLOOD_KERNEL_WEIGHTED,
                           flood_waterlc_and_strm_cells: bool = False,
                           land_cover_water_value: int = 80,
                           find_banks_from_land_cover: bool = False,
                           use_specified_depths: bool = False,
                           specified_depth: float = None,
                           bathy_use_banks: bool = False,
                           disable_bathymetry: bool = False,
                           bathy_args: dict = {},
                           arc_bathymetry_name: str = None,
                           dem_with_bathymetry_name: str = None,
                           overwrite: bool = False) -> Self:
        if not self.reanalysis_flow_file:
            raise ValueError("Base/max flow file must be generated before defining ARC configs.")
        if not self.land_cover:
            raise ValueError("Land cover must be generated before defining ARC configs.")

        self.vdt = os.path.join(self.folder_structure.VDT_Folder, vdt_name)
        self.arc_config = os.path.join(self.folder_structure.ARC_Folder, config_name)
        if curve_file_name:
            self.curve_file = os.path.join(self.folder_structure.VDT_Folder, curve_file_name)
        if ap_file_name:
            self.ap_file = os.path.join(self.folder_structure.VDT_Folder, ap_file_name)
        
        if not disable_bathymetry:
            self.arc_bathymetry = os.path.join(self.folder_structure.bathy_file_folder, arc_bathymetry_name)
            self.dem_with_bathymetry = os.path.join(self.folder_structure.bathy_file_folder, dem_with_bathymetry_name)
        
        if os.path.exists(self.arc_config) and not overwrite:
            return self
        
        os.makedirs(os.path.dirname(self.arc_config), exist_ok=True)
        os.makedirs(os.path.dirname(self.vdt), exist_ok=True)

        params = self._arc_inputs(
            dem=self.cleaned_dem if self.cleaned_dem else self.dem,
            stream_raster=self.stream_raster,
            land_cover=self.land_cover,
            mannings_n_file=self.folder_structure.mannings_text_file,
            reanalysis_flow_file=self.reanalysis_flow_file,
            river_id="COMID",
            highflow_field=maxflow,
            baseflow_field=baseflow,
            bathy_args=bathy_args,
            include_baseflow=not disable_bathymetry,
        )
        
        params["#VDT_Output_File_and_CurveFile"] = ""
        params["VDT_Database_NumIterations"] = bathy_args.get("VDT_Database_NumIterations", 30)
        params["Print_VDT_Database"] = self.vdt
        vdt_suffix = os.path.basename(self.vdt).split("_VDT_Database", 1)[-1]
        params["Print_VDT_Test_File"] = self.vdt.replace(f"_VDT_Database{vdt_suffix}", "_VDT_FS_Bathy.csv")
        if curve_file_name:
            params["Print_Curve_File"] = self.curve_file
        if ap_file_name:
            params["Print_AP_Database"] = self.ap_file
        params["Reach_Average_Curve_File"] = reach_average_curve_file

        if not disable_bathymetry:
            if not self.baseflow_file:
                raise ValueError("Baseflow file must be generated before defining ARC configs with bathymetry.")
            
            params["#Mapper Input Data"] = ""
            params["Comid_Flow_File"] = self.baseflow_file
            params["StrmShp_File"] = self.stream_geometry
            params["mapper"] = mapper.value
            params["FS_ADJUST_FLOW_BY_FRACTION"] = bathy_args.get("FS_ADJUST_FLOW_BY_FRACTION", 1.0)
            params["OutFLD"] = os.path.join(self.folder_structure.flood_folder, floodmap_name)
            params["OutSHP"] = os.path.join(self.folder_structure.flood_folder, floodmap_vector_name)
            params["TopWidthDistanceFactor"] = bathy_args.get("TopWidthDistanceFactor", 1.5)
            params["TW_MultFact"] = bathy_args.get("TW_MultFact", 1.5)
            params["TopWidthPlausibleLimit"] = bathy_args.get("TopWidthPlausibleLimit", 2000)
            if not bathy_args.get("Make_Output_GPKG", True):
                params["Make_Output_GPKG"] = False
            if mapper.is_curve2flood_fldpln_mapper():
                params.update(self._fldpln_inputs())

            if use_specified_depths:
                params["FloodSpreader_SpecifyDepth"] = specified_depth
            else:
                if not self.water_mask:
                    raise ValueError("Water mask must be generated before defining ARC configs without specified depths.")
                params["BathyWaterMask"] = self.water_mask

            params["#Bathymetry_Information"] = ""
            params["Bathy_Trap_H"] = bathy_args.get("Bathy_Trap_H", 0.2)
            params["Bathy_Use_Banks"] = bathy_use_banks

            if flood_waterlc_and_strm_cells:
                params["Flood_WaterLC_and_STRM_Cells"] = True
                params["LAND_WaterValue"] = land_cover_water_value
            if find_banks_from_land_cover:
                params["FindBanksBasedOnLandCover"] = True

            if not arc_bathymetry_name or not dem_with_bathymetry_name:
                raise ValueError("arc_bathymetry_name and dem_with_bathymetry_name must be provided when bathymetry is not disabled.")
            os.makedirs(self.folder_structure.bathy_file_folder, exist_ok=True)
            params["BATHY_Out_File"] = self.arc_bathymetry
            params["AROutBATHY"] = self.arc_bathymetry
            params["FSOutBATHY"] = self.dem_with_bathymetry

        self.write_config(self.arc_config, params)
        
        return self
    
    @ignore_if_dead
    def make_flood_flow_file_from_base_max_file(self, 
                                                    name: str, 
                                                    columns: list[str],  
                                                    overwrite: bool = False, 
                                                    add_to_flow_files: bool = True) -> Self:
        if not self.reanalysis_flow_file:
            raise ValueError("Base/max flow file must be generated before generating flood flow file.")

        flood_flow_file = os.path.join(self.folder_structure.FLOW_Folder, name)
        if add_to_flow_files:
            self.flood_flow_files.append(flood_flow_file)
        if os.path.exists(flood_flow_file) and not overwrite:
            return

        df = pd.read_parquet(self.reanalysis_flow_file) if self.reanalysis_flow_file.endswith('.parquet') else pd.read_csv(self.reanalysis_flow_file)
        df = df.reset_index()[columns]
        df.to_parquet(flood_flow_file, index=False, compression='brotli') if flood_flow_file.endswith('.parquet') else df.to_csv(flood_flow_file, index=False)

        return self
    
    @ignore_if_dead
    def make_baseflow_file_from_base_max_file(self, name: str, columns: list[str], overwrite: bool = False) -> Self:
        self.make_flood_flow_file_from_base_max_file(name=name, columns=columns, overwrite=overwrite, add_to_flow_files=False)
        self.baseflow_file = os.path.join(self.folder_structure.FLOW_Folder, name)

    @ignore_if_dead
    def make_flow_files_from_forecast(self, 
                                          flow_file_prefix: str,
                                          streamflow_source: StreamflowSource,
                                          river_id_field: str,
                                          forecastdate: str = None,
                                          forecasthour: str = None,
                                          nwm_api_key: str = None,
                                          configs: NencartaConfig = None,
                                          overwrite: bool = False) -> Self:
        # now lets download the forecast streamflows
        #Forecast flow data from GeoGLOWS
        # parquet_file_from_geoglows = 'v2-model-table.parquet'     #http://geoglows-v2.s3-website-us-west-2.amazonaws.com/#tables/
        if not self.stream_geometry:
            raise ValueError("Stream geometry must be generated before generating flow files from forecast.")
        
        rivids = Vector(self.stream_geometry).to_geopandas(columns=[river_id_field])[river_id_field].astype(int).tolist()
        
        if forecastdate:
            LOG.info(f"Using forensic forecast date: {forecastdate}")
            if streamflow_source.is_nwm():
                LOG.info(f"Using forensic forecast hour: {forecasthour}")
            else:
                forecasthour = None
            
            if configs.streamflow_source.is_nwm():
                flow_file_name = f'{flow_file_prefix}_{str(forecastdate)}_{forecasthour}_{streamflow_source}_forecast.csv'
            else:
                flow_file_name = f'{flow_file_prefix}_{str(forecastdate)}_{streamflow_source}_forecast.csv'
                
            try:
                ForecastFlowFile = os.path.join(self.folder_structure.FLOW_Folder, flow_file_name)
                if not os.path.exists(flow_file_name) or overwrite:
                    ForecastFlows.Process_and_Write_Forecast_Data(forecastdate, forecasthour, rivids, ForecastFlowFile, streamflow_source, nwm_api_key)
            except Exception as e:
                LOG.error('Could not process forensic forecast streamflow download, please check your date or try again later...')
                raise e
        else:
            # cycle through today through 12 days ago to find the most recent day with a forecast
            found = False
            for fd in range(0,13):
                for fh in range(0,24):
                    try:
                        ForecastFlowFile = os.path.join(self.folder_structure.FLOW_Folder, flow_file_name)
                        forecastdate, forecasthour = ForecastFlows.Get_Date_For_Forecast(fd, fh, streamflow_source) 
                        LOG.info(f"Attempting to download forecast for date: {forecastdate} and hour: {forecasthour}...")
                        # we only need the forecast date for GEOGLOWS, for NWM we need the forecast hour as well               
                        if configs.streamflow_source.is_nwm():
                            flow_file_name = f'{flow_file_prefix}_{str(forecastdate)}_{forecasthour}_{streamflow_source}_forecast.csv'
                        else:
                            flow_file_name = f'{flow_file_prefix}_{str(forecastdate)}_{streamflow_source}_forecast.csv'
                                
                        ForecastFlowFile = os.path.join(self.folder_structure.FLOW_Folder, flow_file_name)
                        if not os.path.exists(ForecastFlowFile) or overwrite:
                            ForecastFlows.Process_and_Write_Forecast_Data(forecastdate, forecasthour, rivids, ForecastFlowFile, streamflow_source, nwm_api_key)
                        
                        if configs:
                            configs.forecastdate = forecastdate
                            configs.forecasthour = forecasthour

                        found = True
                        break
                    except Exception as e:
                        LOG.error(f'Could not process forecast, moving back another day.. ({e})')
                if found:
                    break  # break outer

        LOG.info('Forecast data save here: ' + ForecastFlowFile)    
        self.flood_flow_files.append(ForecastFlowFile)
        return self
        
    @ignore_if_dead
    def define_c2f_configs(self, 
                           config_name: str,
                           flow_file: str,
                           floodmap_name: str,
                           flood_lc_and_streams: bool = False,
                           mapper: Mapper = Mapper.CURVE2FLOOD_KERNEL_WEIGHTED,
                           floodmap_args: dict = {},
                           disable_bathymetry: bool = False,
                           floodmap_vector_name: str = None,
                           velocity_name: str = None,
                           wse_name: str = None,
                           depth_name: str = None,
                           land_watervalue: int = 80,
                           overwrite: bool = False) -> Self:
        if not self.vdt:
            raise ValueError("VDT file must be generated before defining C2F configs.")

        c2f_config = os.path.join(self.folder_structure.ARC_Folder, config_name)
        self.floodmap_configs.append(c2f_config)
        floodmap = os.path.join(self.folder_structure.flood_folder, floodmap_name)
        self.flood_maps.append(floodmap)
        os.makedirs(self.folder_structure.flood_folder, exist_ok=True)

        if floodmap_vector_name:
            floodmap_geometry_file = os.path.join(self.folder_structure.flood_folder, floodmap_vector_name)
            self.floodmap_vectors.append(floodmap_geometry_file)

        if velocity_name:
            velocity_file = os.path.join(self.folder_structure.flood_folder, velocity_name)
            self.velocity_maps.append(velocity_file)

        if wse_name:
            wse_file = os.path.join(self.folder_structure.flood_folder, wse_name)
            self.wse_maps.append(wse_file)

        if depth_name:
            depth_file = os.path.join(self.folder_structure.flood_folder, depth_name)
            self.depth_maps.append(depth_file)

        if os.path.exists(c2f_config) and not overwrite:
            return self
        
        params = {'#ARC_Inputs': ''}
        if disable_bathymetry:
            if self.cleaned_dem:
                params["DEM_File"] = self.cleaned_dem
            else:
                params["DEM_File"] = self.dem
        else:
            if not self.dem_with_bathymetry:
                raise ValueError("DEM with bathymetry must be generated before defining C2F configs with bathymetry.")
            params["DEM_File"] = self.dem_with_bathymetry

        params["Stream_File"] = self.stream_raster
        params["LU_Manning_n"] = self.folder_structure.mannings_text_file

        params["#VDT_Output_File_and_CurveFile"] = ""
        
        params["Print_VDT_Database"] = self.vdt
        if self.curve_file:
            params["Print_Curve_File"] = self.curve_file
        
        params["#Mapper Input Data"] = ""
        params["StrmShp_File"] = self.stream_geometry
        params["Comid_Flow_File"] = flow_file
        params["mapper"] = mapper.value
        if mapper.is_curve2flood_fldpln_mapper():
            params.update(self._fldpln_inputs())
        params["FS_ADJUST_FLOW_BY_FRACTION"] = floodmap_args.get("FS_ADJUST_FLOW_BY_FRACTION", 1.0)
        params["TW_MultFact"] = floodmap_args.get("TW_MultFact", 1.5)
        params["TopWidthPlausibleLimit"] = floodmap_args.get("TopWidthPlausibleLimit", 6000)
        if not floodmap_args.get('Make_Output_GPKG', True):
            params["Make_Output_GPKG"] = False

        if flood_lc_and_streams or velocity_name or True:
            if not self.land_cover:
                raise ValueError("Land cover must be generated before defining C2F configs with flood_lc_and_streams or make_velocity.")
            
            params["Flood_WaterLC_and_STRM_Cells"] = flood_lc_and_streams
            params["LU_Raster_SameRes"] = self.land_cover
            params["LAND_WaterValue"] = land_watervalue

        params["OutFLD"] = floodmap
        if floodmap_vector_name:
            params["OutSHP"] = floodmap_geometry_file
        if velocity_name:
            params["OutVEL"] = velocity_file
        if wse_name:
            params["OutWSE"] = wse_file
        if depth_name:
            params["OutDEP"] = depth_file

        self.write_config(c2f_config, params)

        return self

    @ignore_if_dead
    def make_water_mask_for_dem_cleaner(self,
                                        name: str,
                                        watervalue: int = 80,
                                        overwrite: bool = False) -> Self:
        if not self.land_cover:
            raise ValueError("Land cover must be generated before making water mask for DEM cleaner.")
        if not self.stream_raster:
            raise ValueError("Stream raster must be generated before making water mask for DEM cleaner.")
        
        self.water_mask = os.path.join(self.folder_structure.flood_folder, name)
        if os.path.exists(self.water_mask) and not overwrite:
            return self

        land_cover = Raster(self.land_cover).read_array()
        stream_raster = Raster(self.stream_raster)
        streams = stream_raster.read_array()

        # Mark streams in LC with 1, other areas as -9999
        # Mark streams in SN with 1, other areas with 0
        # Combine LC and SN
        # Mark non-stream areas as -9999 in the final flood map

        F = np.where(
            (streams > 0) | (land_cover == watervalue),
            np.int16(1),
            np.int16(-9999),
        ).reshape(streams.shape)

        Raster.write_array_using_reference(
            F,
            stream_raster,
            self.water_mask,
            gdal.GDT_Int16
        )

    @ignore_if_dead
    def make_clean_dem(self, 
                       mapper: Mapper,
                       flow_file_name: str,
                       use_specified_depth_for_bathy_mask: bool,
                       timer = nullcontext,
                       quiet: bool = False,
                       overwrite: bool = False) -> Self:
        if not self.cleaned_dem:
            raise ValueError("Cleaned DEM must be defined in configs before making clean DEM.")
        if not self.curve_file_for_clean_dem:
            raise ValueError("Curve file for clean DEM must be defined in configs before making clean DEM.")
        if not self.water_mask:
            raise ValueError("Water mask must be generated before making clean DEM.")

        # Run the DEM Cleaner Program, if you wanna
        if os.path.exists(self.cleaned_dem) and not overwrite:
            return

        if overwrite:
            dem_name = os.path.basename(self.dem)
            # Ensure we regenerate ARC + DEM cleaner products when overwriting, so the
            # outputs match a fresh legacy run (important for ignore/*_original fixtures).
            for path in (
                self.cleaned_dem,
                self.curve_file_for_clean_dem,
                getattr(self, "vdt_for_clean_dem", None),
            ):
                if path and os.path.exists(path):
                    try:
                        os.remove(path)
                    except OSError:
                        pass

            for intermediate_name in (
                f"Elev_Streams_{dem_name}",
                f"FLOOD_IMPACT_{dem_name}",
                f"SEED_CONNECT_{dem_name}",
            ):
                intermediate_path = os.path.join(self.folder_structure.dem_updated_folder, intermediate_name)
                if os.path.exists(intermediate_path):
                    os.remove(intermediate_path)
        
        if overwrite or not os.path.exists(self.curve_file_for_clean_dem):
            # start time for the simulation
            with timer('arc_initial'):
                arc = Arc(self.clean_dem_config, quiet=quiet)
                arc.run() 
        else:
            LOG.info(f"{self.curve_file_for_clean_dem} exists and we aren't making it again...")

        with timer('initial_flood_for_cleaner'):
            if mapper == Mapper.FLOODSPREADER and use_specified_depth_for_bathy_mask:
                # Resolve the path to floodspreader.py
                script_dir = os.path.dirname(os.path.abspath(__file__))  # Directory of main.py
                floodspreader_path = os.path.join(script_dir, "floodspreader.py")
                subprocess.run(
                    [sys.executable, floodspreader_path, self.clean_dem_config],
                    check=True,
                )
            elif use_specified_depth_for_bathy_mask:
                LOG.info(f"Executing Curve2Flood using {self.clean_dem_config}")
                Curve2Flood_MainFunction(self.clean_dem_config, quiet=quiet)
        
        OutputID = 'COMID'
        Q_Fraction = 0.10
        TopWidthPlausibleLimit = 600
        search_dist_for_min_elev = 10
        search_dist_perp_cells = 10 # this was 40
        FlowFileName = os.path.join(self.folder_structure.FLOW_Folder, flow_file_name)
        self.make_flood_flow_file_from_base_max_file(
            name=flow_file_name,
            columns=[OutputID, 'p_exceed_50'],
            add_to_flow_files=False
        )

        # start time for the simulation
        with timer('dem_cleaner'):
            DEM_Cleaner.DEM_Cleaner_Program(OutputID, 
                                            self.stream_geometry, 
                                            os.path.dirname(self.dem), 
                                            [os.path.basename(self.dem)], 
                                            [self.stream_raster], 
                                            self.folder_structure.dem_updated_folder, 
                                            FlowFileName, 
                                            self.curve_file_for_clean_dem, 
                                            self.water_mask, 
                                            Q_Fraction, 
                                            TopWidthPlausibleLimit, 
                                            search_dist_for_min_elev, 
                                            search_dist_perp_cells)
