import warnings
from pathlib import Path

import numpy as np
from numba import njit
from osgeo import gdal, ogr

from nencarta.logger import LOG
from nencarta.core.raster import Raster
from nencarta.workspace import Workspace
from nencarta.core.configs import NencartaConfig

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

def _clean_stream_raster(stream_ds: gdal.Dataset) -> None:
    # Get stream raster
    shape = (stream_ds.RasterYSize + 2, stream_ds.RasterXSize + 2)
    array = np.zeros(shape, dtype=np.int64)
    array[1:-1, 1:-1] = stream_ds.ReadAsArray()

    array, first_pass, second_pass = _clean_stream_raster_array(array)
    LOG.info(f'First Pass - Removed {first_pass} cells.')
    LOG.info(f'Second Pass - Removed {second_pass} cells.')
    
    # Write the cleaned array to the raster
    stream_ds.WriteArray(array[1:-1, 1:-1])

def make_stream_raster(workspace: Workspace) -> Path:
    configs = workspace.configs
    if workspace.STRM_File_Clean.exists() and not configs.process_stream_network:
        LOG.info(f"{workspace.STRM_File_Clean} already exists and we aren't making it again...")
        return workspace.STRM_File_Clean
    
    if not workspace.DEM_StrmShp.exists() and not configs.raise_errors_if_nothing_in_domain:
        return None
    
    LOG.info("Creating stream raster from stream geometry...")
    workspace.STRM_File_Clean.parent.mkdir(parents=True, exist_ok=True)

    warnings.filterwarnings(
        "once",
        message="Failed to fetch spatial reference on layer.*",
        category=RuntimeWarning,
    )
    
    dem_raster = Raster(workspace.assigned_dem)
    driver: gdal.Driver = gdal.GetDriverByName('GTiff')

    if workspace.STRM_File_Clean.exists():
        try:
            delete_failed = driver.Delete(workspace.STRM_File_Clean) != 0
        except (PermissionError, RuntimeError, OSError):
            delete_failed = True

        if delete_failed:
            raise PermissionError(
                f"Unable to overwrite stream raster because it is locked: {workspace.STRM_File_Clean}"
            )

    source_ds: gdal.Dataset = gdal.OpenEx(workspace.DEM_StrmShp, gdal.OF_VECTOR)
    if source_ds is None:
        raise RuntimeError(f"Unable to open stream geometry: {workspace.DEM_StrmShp}")

    layer: ogr.Layer = source_ds.GetLayer(0)
    layer_name = layer.GetName() if layer is not None else None

    stream_ds: gdal.Dataset = driver.Create(workspace.STRM_File_Clean, dem_raster.ncols, dem_raster.nrows, 1, gdal.GDT_Int32, ['COMPRESS=DEFLATE', 'PREDICTOR=2'])
    stream_ds.SetGeoTransform(dem_raster.geotransform)
    stream_ds.SetProjection(dem_raster.projection)
    options = gdal.RasterizeOptions(
        attribute=configs.stream_id_field,
        layers=[layer_name] if layer_name and workspace.DEM_StrmShp.suffix.lower().endswith(".gpkg") else None,
    )
    gdal.Rasterize(stream_ds, source_ds, options=options)
    source_ds = None

    # Clean the raster
    _clean_stream_raster(stream_ds)
    stream_ds = None

    return workspace.STRM_File_Clean