import math
from pathlib import Path

import numpy as np
from numba import njit
from osgeo import gdal
from scipy.ndimage import distance_transform_edt, label

from nencarta.logger import LOG
from nencarta.core.raster import Raster
from nencarta.workspace import Workspace
from curve2flood import remove_cells_not_connected

# Land cover's water this far above the river's level on the DEM isn't the river (see near_the_rivers_level)
WATER_MAX_RISE_M = 1.0


def leveled_water(dem: np.ndarray, valid: np.ndarray, min_feature_size: int = 5) -> np.ndarray:
    """The DEM's cells at whole multiples of 0.5 m, as FABDEM levels water, in patches of more than
    ``min_feature_size`` cells (fewer is chance)."""
    leveled = valid & (dem % 0.5 == 0)
    labels, _ = label(leveled, structure=np.ones((3, 3), int))
    counts = np.bincount(labels.ravel())
    counts[0] = 0
    return counts[labels] > min_feature_size


def near_the_rivers_level(dem: np.ndarray, streams: np.ndarray, max_rise: float = WATER_MAX_RISE_M,
                          cells: np.ndarray | None = None) -> np.ndarray:
    """
    The cells of ``dem`` (before the burn) no more than ``max_rise`` above the river's level beside them, where the
    DEM has the river's water. ``streams`` is the stream raster: each stream cell's reach id (> 0), or True.

    The river's level beside a cell is the lowest of the stream cells within twice as far as the nearest one and 5
    cells more, each carried up its reach's slope by how far it lies along the river from the nearest one. Not just
    the nearest one's level: a tributary's line can run across a dry bar beside the river. And not the lowest one's
    as it is: down a sloping river the water beside a cell is higher than the river farther down. Land cover's water
    higher than that is ground the DEM shows dry, a bar or a floodplain the imagery caught flooded, and ARC's
    bathymetry, spread across the water mask, would carve it down to the river's bed.

    Only ``cells`` (default: all) are tested; the rest are False.
    """
    ids = np.asarray(streams).astype(np.int64)
    on_stream = ids > 0
    if not on_stream.any():
        return np.zeros(dem.shape, bool)
    dist, nearest = distance_transform_edt(~on_stream, return_indices=True)
    nearest_z = dem[nearest[0], nearest[1]]
    del nearest
    # The level is never above the nearest stream cell's, so only cells within max_rise of that can pass
    candidates = dem <= nearest_z + max_rise
    if cells is not None:
        candidates &= cells
    del nearest_z

    rows, cols = np.nonzero(on_stream)
    slope = _reach_slopes(dem[rows, cols].astype(np.float64), ids[rows, cols])
    block = 16
    n_block_cols = -(-dem.shape[1] // block)
    bucket = (rows // block) * n_block_cols + cols // block
    order = np.argsort(bucket, kind="mergesort")
    rows, cols, bucket, slope = rows[order], cols[order], bucket[order], slope[order]
    count = np.bincount(bucket, minlength=-(-dem.shape[0] // block) * n_block_cols)
    first = np.concatenate(([0], np.cumsum(count)[:-1]))
    return _near_level(dem.astype(np.float64), np.flatnonzero(candidates), dist, rows, cols,
                       dem[rows, cols].astype(np.float64), slope, first, count, block, n_block_cols, max_rise)


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
def _near_level(dem: np.ndarray, candidates: np.ndarray, dist: np.ndarray, rows: np.ndarray, cols: np.ndarray,
                z: np.ndarray, slope: np.ndarray, first: np.ndarray, count: np.ndarray, block: int,
                n_block_cols: int, max_rise: float) -> np.ndarray:
    """near_the_rivers_level for ``candidates``, with the stream cells bucketed in ``block``-cell squares."""
    nrows, ncols = dem.shape
    n_block_rows = (nrows + block - 1) // block
    out = np.zeros(dem.shape, np.bool_)
    for k in range(candidates.size):
        r = candidates[k] // ncols
        c = candidates[k] % ncols
        d = dist[r, c]
        reach = 2.0 * d + 5.0
        lowest_allowed = dem[r, c] - max_rise
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


def water_mask(workspace: Workspace, dem_raster: Raster, streams: np.ndarray) -> np.ndarray:
    """
    The bathymetry water mask for ``streams`` (reach ids, 0 off them) on ``dem_raster`` (before the burn): the
    stream cells and the water beside them at about the river's level (see near_the_rivers_level), in patches
    touching a stream cell. The water is the land cover's, and with use_dem_derived_channel_mask also the water the
    DEM levels (see leveled_water): the channel this DEM has, where the imagery saw the river elsewhere.
    """
    configs = workspace.configs
    dem = dem_raster.read_array()
    water = np.zeros(dem.shape, bool)
    if workspace.LAND_File.exists():
        water = (Raster(workspace.LAND_File).read_array() == configs.land_watervalue).reshape(dem.shape)
    if configs.use_dem_derived_channel_mask:
        valid = np.isfinite(dem) & (dem != 0)
        if dem_raster.nodata_value is not None:
            valid &= dem != dem_raster.nodata_value
        water |= leveled_water(dem, valid)
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
        stream_ras = Raster(workspace.new_stream_raster)
    else:
        stream_ras = Raster(workspace.STRM_File_Clean)
    streams = stream_ras.read_array()

    F = water_mask(workspace, Raster(workspace.assigned_dem), streams).astype(np.uint8, copy=False)

    Raster.write_array_using_reference(
        F,
        stream_ras,
        workspace.bathy_water_mask,
        gdal.GDT_Byte,
    )

    return workspace.bathy_water_mask
