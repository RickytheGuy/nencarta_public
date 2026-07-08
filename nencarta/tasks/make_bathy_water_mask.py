from pathlib import Path

import numpy as np
from osgeo import gdal

from nencarta.logger import LOG
from nencarta.core.raster import Raster
from nencarta.workspace import Workspace

def make_water_mask(workspace: Workspace) -> Path:
    configs = workspace.configs
    if workspace.bathy_water_mask.exists() and not configs.process_stream_network:
        LOG.info(f"{workspace.bathy_water_mask} already exists and we aren't making it again...")
        return workspace.bathy_water_mask
    
    LOG.info("Creating water mask for DEM cleaner...")
    workspace.bathy_water_mask.parent.mkdir(parents=True, exist_ok=True)
    
    land_cover_array = Raster(workspace.LAND_File).read_array()
    stream_ras = Raster(workspace.STRM_File_Clean)
    streams = stream_ras.read_array()

    # Mark streams in LC with 1, other areas as -9999
    # Mark streams in SN with 1, other areas with 0
    # Combine LC and SN
    # Mark non-stream areas as -9999 in the final flood map

    F = ((streams > 0) | (land_cover_array == configs.land_watervalue)).reshape(streams.shape).astype(np.uint8, copy=False)

    Raster.write_array_using_reference(
        F,
        stream_ras,
        workspace.bathy_water_mask,
        gdal.GDT_Byte,
    )

    return workspace.bathy_water_mask