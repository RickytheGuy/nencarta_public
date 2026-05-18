import os
import tempfile
import datetime
import warnings
from abc import ABC, abstractmethod
from typing import Literal
try:
    from typing import Self
except ImportError:
    from typing_extensions import Self

os.environ["AWS_NO_SIGN_REQUEST"] = "YES"
os.environ["AWS_S3_ENDPOINT"] = "s3.amazonaws.com"
os.environ["KMP_WARNINGS"] = "0"

import numpy as np
import pandas as pd
import geopandas as gpd
from osgeo import gdal, ogr
from numba import njit
from shapely.geometry import box
from propcache import cached_property

from nencarta.logger import LOG
from .raster import Raster
from .vector import Vector

from nencarta._constants import ESA_TILES_FILE

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

@njit(cache=True, nogil=True)
def _clean_stream_raster_numba(array: np.ndarray) -> np.ndarray:
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

    return first_pass, second_pass

def clean_stream_raster(stream_ds: gdal.Dataset, num_passes: int = 2) -> None:
    """
    This function comes from Mike Follum's ARC at https://github.com/MikeFHS/automated-rating-curve
    """
    assert num_passes > 0, "num_passes must be greater than 0"
    
    # Get stream raster
    shape = (stream_ds.RasterYSize + 2, stream_ds.RasterXSize + 2)
    array = np.zeros(shape, dtype=np.int32)
    array[1:-1, 1:-1] = stream_ds.ReadAsArray()

    first_pass, second_pass = _clean_stream_raster_numba(array)
    LOG.info(f'First Pass - Removed {first_pass} cells.')
    LOG.info(f'Second Pass - Removed {second_pass} cells.')
    
    # Write the cleaned array to the raster
    stream_ds.WriteArray(array[1:-1, 1:-1])


class Domain():
    def __init__(self, directory: str):
        self.directory = directory
        self.dem = None
        self.source_dems = []
        self.stream_geometry = None
        self.source_stream_geometry = []
        self.stream_attribute = None
        self.land_cover = None
        self.water_mask = None
        self.stream_raster = None
        self.base_max_flow_file = None
        self.arc_config = None
        self.vdt = None
        self.arc_bathy = None
        self.c2f_baseflow_file = None
        self.burned_dem = None
        self.baseflow_floodmap = None
        self.flood_flow_files = set()
        self.flood_maps = set()

        self.dead = False

        os.makedirs(self.directory, exist_ok=True)

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
            self.directory,
            "dems",
            f"{'buffered_' if buffer else ''}{basename}.{'vrt' if vrt else 'tif'}"
        )

        self.dem = output_dem
        self._id = os.path.splitext(os.path.basename(self.dem))[0]

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
                self._id = os.path.splitext(os.path.basename(dem))[0]
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
    
    @ignore_if_dead
    def generate_stream_raster_from_RFS(self, 
                                        stream_geometry: str | list[str], 
                                        attribute: str = 'LINKNO', 
                                        raise_error_if_no_streams: bool = True,
                                        overwrite: bool = False) -> Self:
        if not self.dem:
            raise ValueError("DEM must be assigned before generating stream raster.")
        
        self.stream_geometry = os.path.join(self.directory, "streams", f"streamlines.parquet")
        self.stream_raster = os.path.join(self.directory, "streams", f"streamlines_{self._id}.tif")
        self.stream_attribute = attribute
        os.makedirs(os.path.dirname(self.stream_geometry), exist_ok=True)

        if not os.path.exists(self.stream_geometry) or not os.path.exists(self.stream_raster) or overwrite:
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

            if gdf.empty:
                if raise_error_if_no_streams:
                    raise ValueError("No stream geometries intersect the DEM extent.")
                
                self.dead = True
                return self

            Vector.save_any_geom(gdf, self.stream_geometry, compression='brotli', write_covering_bbox=True)

        if os.path.exists(self.stream_raster) and not overwrite:
            return self
        
        warnings.filterwarnings(
            "once",
            message="Failed to fetch spatial reference on layer.*",
            category=RuntimeWarning,
        )
        
        dem_raster: Raster = self.dem_raster
        stream_ds: gdal.Dataset = gdal.GetDriverByName('GTiff').Create(self.stream_raster, dem_raster.shape[1], dem_raster.shape[0], 1, gdal.GDT_Int32, options=['COMPRESS=DEFLATE', 'PREDICTOR=2'])
        stream_ds.SetGeoTransform(dem_raster.geotransform)
        stream_ds.SetProjection(dem_raster.projection)

        # Rasterize the streams
        temp: gdal.Dataset = ogr.Open(self.stream_geometry)
        layer = temp.GetLayer()
        gdal.RasterizeLayer(stream_ds, 
                            [1], 
                            layer, 
                            options=[f"ATTRIBUTE={attribute}"],)
        temp = None
        stream_ds.FlushCache()

        # Clean the raster
        clean_stream_raster(stream_ds)
        stream_ds = None

        return self
    
    @ignore_if_dead
    def generate_land_cover(self, land_cover_cache: list[str] = None, vrt: bool = False, overwrite: bool = False) -> Self:
        if not self.dem:
            raise ValueError("DEM must be assigned before generating land cover.")
        
        self.land_cover = os.path.join(self.directory, "land_cover", f"lc_{self._id}.{'vrt' if vrt else 'tif'}")
        if os.path.exists(self.land_cover) and not overwrite:
            return self
        
        os.makedirs(os.path.dirname(self.land_cover), exist_ok=True)

        bounds = self.dem_raster.epsg_4326_bbox
        bbox = box(*bounds)
        tiles = set(gpd.read_file(ESA_TILES_FILE, columns=['ll_tile'], bbox=bbox, ignore_geometry=True, use_arrow=True)['ll_tile'])

        landcover_files = []
        if land_cover_cache:
            cached_files = {os.path.splitext(os.path.basename(path))[0]: path for path in land_cover_cache}
        for tile in tiles:
            if tile in cached_files:
                landcover_files.append(cached_files[tile])
            else:
                landcover_files.append(f"/vsis3/esa-worldcover/v200/2021/map/ESA_WorldCover_10m_2021_v200_{tile}_Map.tif")

        if not landcover_files:
            # Let's make a fake landcover file for arc to use. 
            # Fill it with 10, since the areas that don't have it tend to be tropical (10 is trees)
            ds: gdal.Dataset = gdal.GetDriverByName('GTiff').Create(self.land_cover, self.dem_raster.shape[1], self.dem_raster.shape[0], 1, gdal.GDT_Byte, {'COMPRESS': 'DEFLATE', 'PREDICTOR': '2'})
            ds.SetGeoTransform(self.dem_raster.geotransform)
            ds.SetProjection(self.dem_raster.projection)
            ds.GetRasterBand(1).Fill(10)
            return
        
        proj = self.dem_raster.projection
        if vrt:
            gt = self.dem_raster.geotransform
            xres = gt[1]
            yres = abs(gt[5])
            options = gdal.BuildVRTOptions(outputBounds=bounds,
                                outputSRS=proj,
                                xRes=xres,
                                yRes=yres,
                                resampleAlg='mode')
            gdal.BuildVRT(self.land_cover, landcover_files, options=options)
        else:
            options = gdal.WarpOptions(format='GTiff',
                                outputType=gdal.GDT_Byte,
                                creationOptions=["COMPRESS=ZSTD", "PREDICTOR=2"],
                                outputBounds=bounds,
                                width=self.dem_raster.shape[1],
                                height=self.dem_raster.shape[0],
                                outputBoundsSRS=proj,
                                resampleAlg='mode',
                                dstSRS=proj)
            gdal.Warp(self.land_cover, landcover_files, options=options)

        return self
    
    @ignore_if_dead
    def generate_bathy_water_mask(self, water_class: int = 80, overwrite: bool = False) -> Self:
        if not self.land_cover:
            raise ValueError("Land cover must be generated before generating bathymetry water mask.")
        
        self.water_mask = os.path.join(self.directory, "land_cover", f"water_mask_{self._id}.tif")
        if os.path.exists(self.water_mask) and not overwrite:
            return self
        
        os.makedirs(os.path.dirname(self.water_mask), exist_ok=True)

        lc_ds: gdal.Dataset = gdal.Open(self.land_cover)
        array: np.ndarray = lc_ds.ReadAsArray()
        array = (array == water_class)

        driver = gdal.GetDriverByName('GTiff')
        out_ds: gdal.Dataset = driver.Create(self.water_mask, lc_ds.RasterXSize, lc_ds.RasterYSize, 1, gdal.GDT_Byte, options=['COMPRESS=ZSTD', 'PREDICTOR=2'])
        out_ds.SetGeoTransform(lc_ds.GetGeoTransform())
        out_ds.SetProjection(lc_ds.GetProjection())
        out_ds.WriteArray(array)
        out_ds.FlushCache()
        out_ds = None
        lc_ds = None

        return self

    @ignore_if_dead
    def generate_base_max_flows(self, 
                                parquet: bool = True,
                                overwrite: bool = False) -> Self:
        """
        This function generates a CSV file containing base and maximum flow values for each stream segment in the domain, based on the stream geometry and precomputed flow datasets. The flow values are derived from both the Flow Duration Curve (FDC) and Return Period (RP) datasets, which are accessed via Dask arrays for efficient computation. The resulting CSV file includes columns for various return periods and exceedance probabilities, as well as "premium" flow values calculated as 1.5 times the base flow plus 50.
        This is inspired by nencarta's equivalent function.
        """
        if not self.stream_geometry:
            raise ValueError("Stream geometry must be generated before generating base/max flows.")
        
        self.base_max_flow_file = os.path.join(self.directory, "flow_files", f"base_max.{'parquet' if parquet else 'csv'}")
        if os.path.exists(self.base_max_flow_file) and not overwrite:
            return self

        os.makedirs(os.path.dirname(self.base_max_flow_file), exist_ok=True)

        river_ids = pd.read_parquet(self.stream_geometry, columns=[self.stream_attribute])[self.stream_attribute].unique()
        rp_ds = _get_rp().sel(river_id=river_ids)[['logpearson3', 'gumbel']]
        
        # Convert Xarray to Dask DataFrame and pivot
        rp_df = rp_ds.to_dataframe().reset_index()

        # find the maximum between the gumbel and logpearson3 return periods and label this new column 'return_period_flow'
        rp_df['return_period_flow'] = rp_df[['gumbel', 'logpearson3']].mean(axis=1)

        # keep just the column 'return_period_flow'
        rp_df = rp_df[['river_id', 'return_period', 'return_period_flow']]

        # Convert 'return_period' to category dtype
        rp_df['return_period'] = rp_df['return_period'].astype('category')
        
        # Pivot the table
        rp_pivot_df = rp_df.pivot_table(index='river_id', columns='return_period', values='return_period_flow', aggfunc='mean', observed=False)

        # Rename columns to indicate return periods
        rp_pivot_df = rp_pivot_df.rename(columns={col: f'rp{int(col)}' for col in rp_pivot_df.columns})

        p_exceedance = np.arange(0, 105, 5, dtype=float)
        fdc_ds = _get_fdc().sel(p_exceed=p_exceedance, river_id=river_ids)

        # Convert Xarray to Dask DataFrame
        fdc_df = fdc_ds.to_dataframe().reset_index()

        fdc_pivot = fdc_df.pivot_table(
            index='river_id',
            columns='p_exceed',
            values='hourly_annual',
            aggfunc='mean'
        )
        fdc_pivot = fdc_pivot.rename(columns={p: f"p_exceed_{p}" for p in fdc_pivot.columns})

        final_df = pd.concat([fdc_pivot, rp_pivot_df], axis=1)

        for col in ['p_exceed_0', 'rp100']:
            final_df[f'{col}_premium'] = final_df[col] * 1.5 + 50

        final_df = final_df.reset_index()

        if parquet:
            final_df.round(3).to_parquet(self.base_max_flow_file, compression='brotli', index=False)
        else:
            final_df.round(3).to_csv(self.base_max_flow_file, index=False)

        return self
        
    @ignore_if_dead
    def define_arc_configs(self, 
                           mannings_n_file: str,
                           baseflow: str = 'p_exceed_50',
                           maxflow: str = 'rp100_premium',
                           cross_section_distance: int = 5000,
                           cross_section_wiggle_angle: float = 6.1,
                           cross_section_wiggle_step: float = 1.5,
                           low_spot_range: int = 2,
                           stream_direction_distance: int = 10,
                           stream_slope_distance: int = 10,
                           stream_slope_method: Literal['local_average', 'local_average_corrected', 'end_points'] = 'local_average',
                           vdt_iterations: int = 15,
                           burn_bathymetry: bool = False,
                           bathy_use_banks: bool = False,
                           bathy_banks_from_lc: bool = False,
                           bathy_trap_h: float = 0.2,
                           bathy_flow_multiplier: float = 1.0,
                           bathy_topwidth_multiplier: float = 1.0,
                           bathy_topwidth_limit: float = 1000,
                           overwrite: bool = False) -> Self:
        if not self.base_max_flow_file:
            raise ValueError("Base/max flow file must be generated before defining ARC configs.")
        if not self.land_cover:
            raise ValueError("Land cover must be generated before defining ARC configs.")

        self.vdt = os.path.join(self.directory, "vdt", f"vdt_{self._id}.parquet")
        self.arc_config = os.path.join(self.directory, "configs", f"arc_{self._id}.tsv")
        self.c2f_baseflow_file = os.path.join(self.directory, "flow_files", f"baseflow.csv")

        if burn_bathymetry and (not os.path.exists(self.c2f_baseflow_file) or overwrite):
            baseflow_df = pd.read_parquet(self.base_max_flow_file) if self.base_max_flow_file.endswith('.parquet') else pd.read_csv(self.base_max_flow_file)
            baseflow_df = baseflow_df[['river_id', baseflow]]
            baseflow_df.to_csv(self.c2f_baseflow_file, index=False)

        if os.path.exists(self.arc_config) and not overwrite:
            return self
        
        os.makedirs(os.path.dirname(self.arc_config), exist_ok=True)
        os.makedirs(os.path.dirname(self.vdt), exist_ok=True)
        
        arc_args = {
            "DEM_File": self.dem,
            "Stream_File": self.stream_raster,
            "LU_Raster_SameRes": self.land_cover,
            "LU_Manning_n": mannings_n_file,
            "Flow_File": self.base_max_flow_file,
            "Flow_File_ID": "river_id",
            "Flow_File_BF": baseflow,
            "Flow_File_QMax": maxflow,
            "X_Section_Dist": cross_section_distance,
            "Degree_Manip": cross_section_wiggle_angle,
            "Degree_Interval":	cross_section_wiggle_step,
            "Low_Spot_Range":  low_spot_range,
            "Gen_Dir_Dist":	stream_direction_distance,
            "Gen_Slope_Dist":	stream_slope_distance,
            "Stream_Slope_Method":	stream_slope_method,
            "VDT_Database_NumIterations":	vdt_iterations, 
            "Print_VDT_Database": self.vdt,   
        }

        if stream_slope_method == "end_points":
            arc_args["StrmShp_File"] = self.stream_geometry
        
        if burn_bathymetry:
            if not self.water_mask:
                raise ValueError("Water mask must be generated before defining ARC configs with bathymetry.")
            
            self.arc_bathy = os.path.join(self.directory, "bathymetry", f"arcbathy_{self._id}.tif")
            self.burned_dem = os.path.join(self.directory, "bathymetry", f"burned_{self._id}.tif")
            self.baseflow_floodmap = os.path.join(self.directory, "flood_maps", f"{baseflow}_{self._id}.tif")
            os.makedirs(os.path.dirname(self.arc_bathy), exist_ok=True)
            os.makedirs(os.path.dirname(self.burned_dem), exist_ok=True)
            os.makedirs(os.path.dirname(self.baseflow_floodmap), exist_ok=True)
            
            arc_args.update({
                "BathyWaterMask":	self.water_mask,
                "Bathy_Trap_H": bathy_trap_h,
                "Bathy_Use_Banks":	bathy_use_banks,
                "FindBanksBasedOnLandCover": bathy_banks_from_lc,
                "AROutBATHY":	self.arc_bathy,
                "BATHY_Out_File":	self.arc_bathy,
                "FSOutBATHY":	self.burned_dem,
                "Comid_Flow_File":	self.c2f_baseflow_file,
                "FS_ADJUST_FLOW_BY_FRACTION":	bathy_flow_multiplier,
                "OutFLD":	self.baseflow_floodmap,
                "TW_MultFact":	bathy_topwidth_multiplier,
                "TopWidthPlausibleLimit":	bathy_topwidth_limit,
                "Make_Output_GPKG":	False,
            })

        with open(self.arc_config, 'w') as f:
            f.write(f"# ARC Config File for {self._id} ({datetime.datetime.now()})\n\n")
            for key, value in arc_args.items():
                f.write(f"{key}\t{value}\n")
        
        return self
    
    @ignore_if_dead
    def generate_flood_flow_file_from_base_max_file(self, columns: str | list[str], parquet: bool = True, overwrite: bool = False) -> Self:
        if not self.base_max_flow_file:
            raise ValueError("Base/max flow file must be generated before generating flood flow file.")

        if isinstance(columns, str):
            columns = [columns]

        flood_flow_file = os.path.join(self.directory, "flow_files", f"{'_'.join(columns)}.{'parquet' if parquet else 'csv'}")
        if os.path.exists(flood_flow_file) and not overwrite:
            return

        df = pd.read_parquet(self.base_max_flow_file) if self.base_max_flow_file.endswith('.parquet') else pd.read_csv(self.base_max_flow_file)
        df = df.reset_index()[['river_id', *columns]]
        df.to_parquet(flood_flow_file, index=False, compression='brotli') if flood_flow_file.endswith('.parquet') else df.to_csv(flood_flow_file, index=False)

        self.flood_flow_files.add(flood_flow_file)
        return self
    
    @ignore_if_dead
    def define_c2f_configs(self, 
                           use_burned_dem: bool = False,
                           flood_lc_and_streams: bool = False,
                           flood_local: bool = False,
                           mapper: Literal['Curve2Flood-Kernel Weighted', 'Curve2Flood-Multi-Point Interpolation'] = 'Curve2Flood-Kernel Weighted',
                           flow_multiplier: float = 1.0,
                           topwidth_multiplier: float = 1.0,
                           topwidth_limit: float = 1000,
                           overwrite: bool = False) -> Self:
        if not self.vdt:
            raise ValueError("VDT file must be generated before defining C2F configs.")
        if use_burned_dem and not self.burned_dem:
            raise ValueError("Burned DEM must be generated before defining C2F configs with burned DEM.")

        for flow_file in self.flood_flow_files:
            c2f_config = os.path.join(self.directory, "configs", f"c2f_{self._id}_{os.path.splitext(os.path.basename(flow_file))[0]}.tsv")
            floodmap = os.path.join(self.directory, "flood_maps", f"{os.path.splitext(os.path.basename(flow_file))[0]}_{self._id}.tif")
            self.flood_maps.add(floodmap)
            if os.path.exists(c2f_config) and not overwrite:
                continue
            
            os.makedirs(os.path.dirname(c2f_config), exist_ok=True)

            c2f_args = {
                "DEM_File":	self.burned_dem if use_burned_dem else self.dem,
                "Stream_File":	self.stream_raster,
                "Print_VDT_Database":	self.vdt,
                "Comid_Flow_File":	flow_file,
                "FS_ADJUST_FLOW_BY_FRACTION":	flow_multiplier,
                "TW_MultFact":	topwidth_multiplier,
                "TopWidthPlausibleLimit":	topwidth_limit,
                "Make_Output_GPKG":	False,
                "Flood_WaterLC_and_STRM_Cells":	flood_lc_and_streams,
                "LU_Raster_SameRes": self.land_cover,
                "LAND_WaterValue":	80,
                "OutFLD":	floodmap,
                "LocalFloodOption":	flood_local,
                "mapper": mapper
            }
            with open(c2f_config, 'w') as f:
                f.write(f"# Curve2Flood Config File for {self._id} ({datetime.datetime.now()})\n\n")
                for key, value in c2f_args.items():
                    f.write(f"{key}\t{value}\n")

        return self
