import math
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from numba import njit
from osgeo import gdal
from scipy.ndimage import distance_transform_edt

from nencarta.logger import LOG
from nencarta.core.raster import Raster
from nencarta.core.vector import Vector
from nencarta.tasks.channel_detection import _CHUNK_ROWS, _feature_distance, river_channel
from nencarta.tasks.stream_conflation import cell_sizes_m
from nencarta.workspace import Workspace
from curve2flood import remove_cells_not_connected

# Land cover's water this far above the river's level on the DEM isn't the river (see near_the_rivers_level)
WATER_MAX_RISE_M = 1.0


def reach_areas_km2(streams_gdf: gpd.GeoDataFrame | pd.DataFrame, configs) -> dict:
    """Each reach's drainage area in km2, by id, from area_km2_field or else area_m2_field (empty without either)."""
    id_col = configs.stream_id_field
    if streams_gdf is None or id_col not in streams_gdf.columns:
        return {}
    if configs.area_km2_field and configs.area_km2_field in streams_gdf.columns:
        area = pd.to_numeric(streams_gdf[configs.area_km2_field], errors="coerce")
    elif configs.area_m2_field and configs.area_m2_field in streams_gdf.columns:
        area = pd.to_numeric(streams_gdf[configs.area_m2_field], errors="coerce") / 1e6
    else:
        return {}
    ids = pd.to_numeric(streams_gdf[id_col], errors="coerce")
    keep = ids.notna() & area.notna()
    return dict(zip(ids[keep].astype(np.int64).tolist(), area[keep].astype(float).tolist()))


def dem_channel(dem_raster, dem: np.ndarray, valid: np.ndarray, streams: np.ndarray, area_km2: dict,
                connected: bool = True, flat_water: bool = True, channel: bool = True):
    """
    The channel the DEM shows along ``streams`` (see channel_detection.river_channel, which skips ``flat_water`` or
    ``channel`` when told to), on the grid of ``dem_raster`` (a Raster, or a gdal.Dataset).
    """
    if isinstance(dem_raster, gdal.Dataset):
        geotransform, projection = dem_raster.GetGeoTransform(), dem_raster.GetProjection()
    else:
        geotransform, projection = dem_raster.geotransform, dem_raster.projection
    dx, dy = cell_sizes_m(geotransform, projection, dem.shape[0])
    return river_channel(dem, valid, streams, area_km2, float(np.median(dx)), dy, connected=connected,
                         flat_water=flat_water, channel=channel)


def near_the_rivers_level(dem: np.ndarray, streams: np.ndarray, max_rise: float = WATER_MAX_RISE_M,
                          cells: np.ndarray | None = None) -> np.ndarray:
    """
    The cells of ``dem`` (before the burn) no more than ``max_rise`` above the river's level beside them, where the
    DEM has the river's water. ``streams`` is the stream raster: each stream cell's reach id (> 0), or True.

    The river's level beside a cell is the lowest of the stream cells within twice as far as the nearest one and 5
    cells more, each carried up its reach's slope by how far it lies along the river from the nearest one. Not just
    the nearest one's level: a tributary's line can run across a dry bar beside the river. And not the lowest one's
    as it is: down a sloping river the water beside a cell is higher than the river farther down. Land cover's water
    higher than that is ground the DEM shows dry, a bar or a floodplain the imagery caught flooded, and ARC, whose
    banks are where the water mask's water ends, would carve it down to the river's bed.

    Only ``cells`` (default: all) are tested; the rest are False.
    """
    streams = np.asarray(streams)
    on_stream = streams > 0
    if not on_stream.any():
        return np.zeros(dem.shape, bool)
    nearest = distance_transform_edt(~on_stream, return_distances=False, return_indices=True)
    # The level is never above the nearest stream cell's, so only cells within max_rise of that can pass. Found a
    # block of rows at a time, with their distance to the nearest stream cell, to keep the temporary arrays small.
    ncols = dem.shape[1]
    found, distances = [], []
    for r0 in range(0, dem.shape[0], _CHUNK_ROWS):
        r1 = min(r0 + _CHUNK_ROWS, dem.shape[0])
        rows, cols = np.nonzero(np.ones((r1 - r0, ncols), bool) if cells is None else cells[r0:r1])
        rows += r0
        near_rows, near_cols = nearest[0][rows, cols], nearest[1][rows, cols]
        keep = dem[rows, cols] <= dem[near_rows, near_cols] + max_rise
        rows, cols = rows[keep], cols[keep]
        found.append(rows * ncols + cols)
        distances.append(_feature_distance(near_rows[keep] - rows, near_cols[keep] - cols, (1.0, 1.0)))
    del nearest
    candidates, distance = np.concatenate(found), np.concatenate(distances)

    rows, cols = np.nonzero(on_stream)
    slope = _reach_slopes(dem[rows, cols].astype(np.float64), streams[rows, cols].astype(np.int64))
    block = 16
    n_block_cols = -(-dem.shape[1] // block)
    bucket = (rows // block) * n_block_cols + cols // block
    order = np.argsort(bucket, kind="mergesort")
    rows, cols, bucket, slope = rows[order], cols[order], bucket[order], slope[order]
    count = np.bincount(bucket, minlength=-(-dem.shape[0] // block) * n_block_cols)
    first = np.concatenate(([0], np.cumsum(count)[:-1]))
    return _near_level(dem.shape[0], ncols, candidates, distance, dem.ravel()[candidates].astype(np.float64), rows,
                       cols, dem[rows, cols].astype(np.float64), slope, first, count, block, n_block_cols, max_rise)


def _reach_slopes(z: np.ndarray, ids: np.ndarray) -> np.ndarray:
    """
    For each stream cell (height ``z``, reach ``ids``), how much its reach falls per cell along it: the spread
    between the 10th and 90th percentiles of its cells' heights over its number of cells.
    """
    order = np.argsort(ids, kind="mergesort")
    bounds = np.flatnonzero(np.diff(ids[order])) + 1
    slope = np.empty(z.size, np.float64)
    for part in np.split(order, bounds):
        heights = z[part]
        slope[part] = max(float(np.percentile(heights, 90) - np.percentile(heights, 10)), 0.0) / part.size
    return slope


@njit(cache=True)
def _near_level(nrows: int, ncols: int, candidates: np.ndarray, distance: np.ndarray, height: np.ndarray,
                rows: np.ndarray, cols: np.ndarray, z: np.ndarray, slope: np.ndarray, first: np.ndarray,
                count: np.ndarray, block: int, n_block_cols: int, max_rise: float) -> np.ndarray:
    """
    near_the_rivers_level for ``candidates`` (flat indices on an ``nrows`` by ``ncols`` grid, each ``distance`` cells
    from the nearest stream cell and ``height`` high), with the stream cells bucketed in ``block``-cell squares.
    """
    n_block_rows = (nrows + block - 1) // block
    out = np.zeros((nrows, ncols), np.bool_)
    for k in range(candidates.size):
        r = candidates[k] // ncols
        c = candidates[k] % ncols
        d = distance[k]
        reach = 2.0 * d + 5.0
        lowest_allowed = height[k] - max_rise
        near = True
        for br in range(max(int((r - reach) // block), 0), min(int((r + reach) // block), n_block_rows - 1) + 1):
            for bc in range(max(int((c - reach) // block), 0), min(int((c + reach) // block), n_block_cols - 1) + 1):
                b = br * n_block_cols + bc
                for i in range(first[b], first[b] + count[b]):
                    d2 = (rows[i] - r) ** 2 + (cols[i] - c) ** 2
                    if d2 > reach * reach:
                        continue
                    if z[i] + slope[i] * math.sqrt(max(d2 - d * d, 0.0)) < lowest_allowed:
                        near = False
                        break
                if not near:
                    break
            if not near:
                break
        out[r, c] = near
    return out


def water_mask(workspace: Workspace, dem_raster: Raster, streams: np.ndarray, area_km2: dict | None = None
               ) -> np.ndarray:
    """
    The bathymetry water mask for ``streams`` (reach ids, 0 off them, draining ``area_km2`` by id) on ``dem_raster``
    (before the burn): the stream cells and the water beside them at about the river's level (see
    near_the_rivers_level), in patches touching a stream cell. The water is the land cover's, and with
    use_dem_derived_channel_mask also the channel the DEM shows (see channel_detection), where the imagery saw the
    river elsewhere.
    """
    configs = workspace.configs
    dem = dem_raster.read_array()
    water = np.zeros(dem.shape, bool)
    if workspace.LAND_File.exists():
        water = (Raster(workspace.LAND_File).read_array() == configs.land_watervalue).reshape(dem.shape)
    if configs.use_dem_derived_channel_mask:
        valid = np.isfinite(dem)
        if dem_raster.nodata_value is not None:
            valid &= dem != dem_raster.nodata_value
        water |= dem_channel(dem_raster, dem, valid, streams, area_km2 or {}, flat_water=False).channel
    water &= near_the_rivers_level(dem, streams, cells=water)
    return remove_cells_not_connected(water | (streams > 0), streams) > 0


def make_water_mask(workspace: Workspace) -> Path:
    """
    Write the water mask (workspace.bathy_water_mask): 1 on stream cells and the water beside them at the river's
    level (see water_mask), 0 elsewhere.
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
    
    # ARC carves the bathymetry along the streams it is given, which are the moved ones when the network has been
    # moved to the thalweg, and curve2flood discards whatever bathymetry falls outside this mask. On an overwrite
    # they are moved again, so the last run's are stale; the mask is remade from the new ones once they are.
    if configs.move_stream_network_to_thalweg and workspace.new_stream_raster.exists() and not configs.overwrite:
        stream_ras, stream_lines = Raster(workspace.new_stream_raster), workspace.new_StrmShp_matched
    else:
        stream_ras, stream_lines = Raster(workspace.STRM_File_Clean), workspace.DEM_StrmShp
    streams = stream_ras.read_array()
    area_km2 = {}
    if configs.use_dem_derived_channel_mask and stream_lines.exists():
        table = Vector(stream_lines, use_threads=not configs.parallel).to_geopandas()
        area_km2 = reach_areas_km2(table, configs)

    F = water_mask(workspace, Raster(workspace.assigned_dem), streams, area_km2).astype(np.uint8, copy=False)

    Raster.write_array_using_reference(
        F,
        stream_ras,
        workspace.bathy_water_mask,
        gdal.GDT_Byte,
    )

    return workspace.bathy_water_mask
