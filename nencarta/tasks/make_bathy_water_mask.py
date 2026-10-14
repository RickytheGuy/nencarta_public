from pathlib import Path

import numpy as np
from osgeo import gdal
from scipy.ndimage import distance_transform_edt, minimum_filter

from nencarta.logger import LOG
from nencarta.core.raster import Raster
from nencarta.workspace import Workspace

# Land cover's water this far above the river's level on the DEM isn't the river (see near_the_rivers_level)
WATER_MAX_RISE_M = 1.0


def near_the_rivers_level(dem: np.ndarray, streams: np.ndarray, max_rise: float = WATER_MAX_RISE_M) -> np.ndarray:
    """
    The cells of ``dem`` (before the burn) no more than ``max_rise`` above the river's level beside them, where the
    DEM has the river's water: the lowest of the ``streams`` cells within twice as far as the nearest one and 5
    cells more (in square windows of up to 256 cells; the nearest one's, farther than that). Not just the nearest
    one's: a tributary's line can run across a dry bar beside the river. Land cover's water higher than that is
    ground the DEM shows dry, a bar or a floodplain the imagery caught flooded, and ARC's bathymetry, spread across
    the water mask, would carve it down to the river's bed.
    """
    if not streams.any():
        return np.zeros(dem.shape, bool)
    distance, nearest = distance_transform_edt(~streams, return_indices=True)
    level = dem[nearest[0], nearest[1]].astype(np.float64)
    del nearest
    on_streams = np.where(streams, dem, np.inf)
    reach = 2 * distance + 5
    del distance
    done = np.zeros(dem.shape, bool)
    for half in (8, 16, 32, 64, 128, 256):
        use = ~done & (reach <= half)
        if use.any():
            level[use] = np.minimum(level[use], minimum_filter(on_streams, size=2 * half + 1)[use])
            done |= use
    return dem <= level + max_rise


def make_water_mask(workspace: Workspace) -> Path:
    """
    Write the water mask (workspace.bathy_water_mask): 1 on stream cells and on water land cover near the river's
    level (see near_the_rivers_level), 0 elsewhere.
    """
    configs = workspace.configs
    if workspace.bathy_water_mask.exists() and not configs.overwrite:
        LOG.info(f"{workspace.bathy_water_mask} already exists and we aren't making it again...")
        return workspace.bathy_water_mask

    # make_stream_raster() makes no stream raster when there are no streams in the domain
    if not workspace.DEM_StrmShp.exists() and not configs.raise_errors_if_nothing_in_domain:
        return None

    LOG.info("Creating water mask...")
    workspace.bathy_water_mask.parent.mkdir(parents=True, exist_ok=True)
    
    land_cover_array = Raster(workspace.LAND_File).read_array()
    # ARC carves the bathymetry along the streams it is given, which are the moved ones when the network has been
    # moved to the thalweg, and curve2flood discards whatever bathymetry falls outside this mask. On an overwrite
    # they are moved again, so the last run's are stale; the mask is remade from the new ones once they are.
    if configs.move_stream_network_to_thalweg and workspace.new_stream_raster.exists() and not configs.overwrite:
        stream_ras = Raster(workspace.new_stream_raster)
    else:
        stream_ras = Raster(workspace.STRM_File_Clean)
    streams = stream_ras.read_array()

    water = (land_cover_array == configs.land_watervalue).reshape(streams.shape)
    water &= near_the_rivers_level(Raster(workspace.assigned_dem).read_array(), streams > 0)
    F = ((streams > 0) | water).astype(np.uint8, copy=False)

    Raster.write_array_using_reference(
        F,
        stream_ras,
        workspace.bathy_water_mask,
        gdal.GDT_Byte,
    )

    return workspace.bathy_water_mask
