"""
Burning the source streams into the DEM, and moving them onto its D8 flow paths (see stream_conflation).
"""
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
from numba import njit
import geopandas as gpd
from osgeo import gdal, ogr, osr
import richdem as rd
from scipy.ndimage import binary_dilation, label, maximum_filter, minimum_filter
from shapely.geometry import Point, LineString

from nencarta.logger import LOG
from nencarta.core.raster import Raster
from nencarta.core.vector import Vector
from nencarta.tasks.make_bathy_water_mask import dem_channel, reach_areas_km2, water_mask
from nencarta.tasks.make_stream_geometry import _filter_streams_by_stream_order
from nencarta.tasks.stream_conflation import (
    D8_OFFSETS,
    ConflationSettings,
    cell_sizes_m,
    conflate_to_flow_directions,
    d8_offset_tables as _d8_offset_tables,
)
from nencarta.tasks.stream_relocation import breach_dammed_rivers, relocate_misplaced_reaches
from nencarta.workspace import Workspace

def _in_memory_ogr_layer(gdf: gpd.GeoDataFrame, projection: str,
                         values: pd.Series | None = None) -> tuple[ogr.DataSource, ogr.Layer]:
    """
    Copy ``gdf``'s geometries, in order, into an in-memory OGR layer for GDAL to rasterize, with ``values`` in an
    integer field named "value" if given. Keep the returned DataSource alive while using the layer.
    """
    # 'Memory' is the pre-GDAL-3.11 spelling and now warns on every call; 'MEM' is the
    # current name. Fall back so this still works on older GDAL builds.
    driver = ogr.GetDriverByName('MEM') or ogr.GetDriverByName('Memory')
    ogr_ds: ogr.DataSource = driver.CreateDataSource('lakes')
    srs = osr.SpatialReference()
    srs.ImportFromWkt(projection)
    layer: ogr.Layer = ogr_ds.CreateLayer('lakes', srs, ogr.wkbUnknown)
    if values is not None:
        layer.CreateField(ogr.FieldDefn("value", ogr.OFTInteger64))
        values = pd.to_numeric(values, errors="coerce").to_numpy()
    layer_defn = layer.GetLayerDefn()
    for i, geometry in enumerate(gdf.geometry.values):
        if geometry is None or geometry.is_empty:
            continue
        feature = ogr.Feature(layer_defn)
        feature.SetGeometry(ogr.CreateGeometryFromWkb(geometry.wkb))
        if values is not None and np.isfinite(values[i]):
            feature.SetField("value", int(values[i]))
        layer.CreateFeature(feature)
        feature = None
    return ogr_ds, layer

def load_lake_array(workspace: Workspace, dem_raster: Raster) -> np.ndarray | None:
    """
    The lakes on the DEM's grid, True in a lake, or None without ``lakes``. Rasterized once to
    workspace.lake_raster. The lakes are read with geopandas, since GDAL may not read GeoParquet.
    """
    if not workspace.configs.lakes:
        return None

    if workspace.lake_raster.exists() and not workspace.configs.overwrite:
        # A boolean mask, as below: curve2flood indexes the DEM with it, and a 0/1 integer array would
        # index rows 0 and 1 instead of the lake cells
        lakes = gdal.Open(str(workspace.lake_raster)).ReadAsArray().astype(np.bool_, copy=False)
        return lakes

    lakes_ds: gdal.Dataset = gdal.GetDriverByName('GTiff').Create(str(workspace.lake_raster), dem_raster.shape[1], dem_raster.shape[0], 1, gdal.GDT_Byte, options=[f'COMPRESS={workspace.configs.compression}'])
    lakes_ds.SetGeoTransform(dem_raster.geotransform)
    lakes_ds.SetProjection(dem_raster.projection)

    # Vector.to_geopandas() does the bbox reprojection itself, so the subset comes back
    # already clipped to the DEM's footprint.
    lakes_gdf = load_lake_gdf(workspace, dem_raster)
    if lakes_gdf is not None and not lakes_gdf.empty:
        lakes_gdf = lakes_gdf.to_crs(dem_raster.projection)
        ogr_ds, lakes_layer = _in_memory_ogr_layer(lakes_gdf, dem_raster.projection)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gdal.RasterizeLayer(lakes_ds, [1], lakes_layer, burn_values=[1])

        ogr_ds = None

    lakes_ds.FlushCache()
    lakes = lakes_ds.ReadAsArray().astype(np.bool_, copy=False)
    return lakes

def load_lake_gdf(workspace: Workspace, dem_raster: Raster) -> gpd.GeoDataFrame | None:
    """The lakes overlapping the DEM, or None without ``lakes``."""
    if not workspace.configs.lakes:
        return None

    lakes_gdf = Vector(workspace.configs.lakes, not workspace.configs.parallel).to_geopandas(bbox_epsg_4326=dem_raster.epsg_4326_bbox)

    return lakes_gdf

# The burn fills the cells between a line's vertices up to this many cells apart (see _burn_linestring)
MAX_BURN_GAP_CELLS = 3

# Flow directions are D8 pointer codes (see stream_conflation.D8_OFFSETS), Int16 with this nodata
FLOWDIR_NODATA = -32768
FILLED_NODATA = -9999.0
_D8_ROWS = np.array([offset[0] for offset in D8_OFFSETS.values()], np.int64)
_D8_COLS = np.array([offset[1] for offset in D8_OFFSETS.values()], np.int64)
_D8_CODES = np.array(list(D8_OFFSETS), np.int16)


def fill_and_flow_directions(dem: np.ndarray, valid: np.ndarray, geotransform: tuple,
                             projection: str) -> tuple[np.ndarray, np.ndarray]:
    """
    ``dem`` with its depressions filled (NaN off ``valid``), and its D8 flow directions (see
    stream_conflation.D8_OFFSETS; 0 where nothing around is lower, FLOWDIR_NODATA off the DEM), with RichDEM.

    Every depression is filled to its spill level (Priority-Flood, Barnes et al. 2014), and the flats that leaves --
    and the DEM's own, as FABDEM's leveled water -- are given a gradient for the flow directions alone (Barnes et al.
    2014, "An efficient assignment of drainage direction over flat surfaces"): away from the higher ground around
    them and toward their outlets, so flow paths run down a flat's middle rather than straight along one edge. Off
    the DEM, its edges and nodata, is lower than anything on it, so a flat running off it drains there. The gradient
    is imposed in double precision: in single precision its increments round away on a large flat, leaving pits.
    """
    valid = np.asarray(valid, bool)
    low = float(np.min(dem[valid])) - 10.0 if valid.any() else 0.0
    padded = np.pad(np.where(valid, dem, low).astype(np.float64), 1, constant_values=low)
    gt = list(geotransform)
    gt[0] -= gt[1]
    gt[3] -= gt[5]
    filled = rd.FillDepressions(rd.rdarray(padded, no_data=-np.finfo(np.float64).max, geotransform=gt),
                                in_place=False)
    routed = np.asarray(rd.ResolveFlats(filled, in_place=False))[1:-1, 1:-1]
    filled = np.asarray(filled)[1:-1, 1:-1]
    dx, dy = cell_sizes_m(geotransform, projection, dem.shape[0])
    flowdir = _d8_pointers(routed, valid, np.asarray(dx, np.float64), float(dy), _D8_ROWS, _D8_COLS, _D8_CODES)
    return np.where(valid, filled, np.nan), flowdir


@njit(cache=True)
def _d8_pointers(z: np.ndarray, valid: np.ndarray, dx_by_row: np.ndarray, dy: float, d_rows: np.ndarray,
                 d_cols: np.ndarray, codes: np.ndarray) -> np.ndarray:
    """The code of each valid cell's steepest descent to a valid neighbour (distances in metres), or 0."""
    nrows, ncols = z.shape
    out = np.full((nrows, ncols), FLOWDIR_NODATA, np.int16)
    for row in range(nrows):
        for col in range(ncols):
            if not valid[row, col]:
                continue
            steepest = 0.0
            code = 0
            for k in range(d_rows.size):
                r, c = row + d_rows[k], col + d_cols[k]
                if r < 0 or r >= nrows or c < 0 or c >= ncols or not valid[r, c]:
                    continue
                slope = (z[row, col] - z[r, c]) / np.sqrt((dx_by_row[row] * d_cols[k]) ** 2 + (dy * d_rows[k]) ** 2)
                if slope > steepest:
                    steepest = slope
                    code = codes[k]
            out[row, col] = code
    return out


def derive_flow_directions(workspace: Workspace, dem_path: Path, filled_path: Path | None = None,
                           flowdir_path: Path | None = None) -> None:
    """
    Fill the depressions in ``dem_path`` and derive its D8 flow directions (see fill_and_flow_directions), writing
    ``filled_path`` and ``flowdir_path`` (default: workspace.filled_dem and workspace.flowdir) with the configs'
    compression.
    """
    raster = Raster(dem_path)
    dem = raster.read_array().astype(np.float64)
    valid = np.isfinite(dem)
    if raster.nodata_value is not None:
        valid &= dem != raster.nodata_value
    filled, flowdir = fill_and_flow_directions(dem, valid, raster.geotransform, raster.projection)
    compression = workspace.configs.compression
    _write_raster(Path(filled_path or workspace.filled_dem), np.where(valid, filled, FILLED_NODATA), raster,
                  gdal.GDT_Float64, FILLED_NODATA, compression)
    _write_raster(Path(flowdir_path or workspace.flowdir), flowdir, raster, gdal.GDT_Int16, FLOWDIR_NODATA,
                  compression)


def _write_raster(path: Path, values: np.ndarray, reference: Raster, data_type: int, nodata: float,
                  compression: str) -> None:
    """Write ``values`` as a tiled GeoTIFF on ``reference``'s grid."""
    path.parent.mkdir(parents=True, exist_ok=True)
    driver: gdal.Driver = gdal.GetDriverByName("GTiff")
    if path.exists():
        driver.Delete(str(path))
    ds: gdal.Dataset = driver.Create(str(path), values.shape[1], values.shape[0], 1, data_type,
                                     [f"COMPRESS={compression}", "TILED=YES", "BIGTIFF=IF_SAFER"])
    ds.SetGeoTransform(reference.geotransform)
    ds.SetProjection(reference.projection)
    band = ds.GetRasterBand(1)
    band.WriteArray(values)
    band.SetNoDataValue(nodata)
    band = None
    # Closed explicitly, as in _write_stream_raster
    ds.Close()

def burn_streams_and_move_streams(workspace: Workspace) -> Path:
    """
    Burn the streams into the DEM (burn_streams), and move them onto its D8 flow paths
    (move_stream_network_to_thalweg). Writes the burned DEM, the filled DEM and flow directions, the moved stream
    network and its raster, and the water mask, plus the stream info table for Curve2Flood-FLDPLNpy.
    """
    if not workspace.DEM_StrmShp.exists() and not workspace.configs.raise_errors_if_nothing_in_domain:
        return None

    configs = workspace.configs
    should_burn_streams = configs.burn_streams and (not workspace.fixed_dem.exists() or configs.overwrite)
    should_move_streams = configs.move_stream_network_to_thalweg and \
        (not workspace.new_StrmShp_matched.exists() or not workspace.new_stream_raster.exists() or configs.overwrite or 
         (configs.mapper.is_curve2flood_fldpln_mapper() and (
            not workspace.stream_info_file.exists() or not workspace.filled_dem.exists() or not workspace.flowdir.exists()
         )))

    assigned_dem = Raster(workspace.assigned_dem)
    flow_directions_ready = False

    if should_burn_streams:
        workspace.dem_updated_folder.mkdir(parents=True, exist_ok=True)

        lakes = load_lake_array(workspace, assigned_dem)

        source_gdf = Vector(workspace.DEM_StrmShp, not configs.parallel).to_geopandas().to_crs(assigned_dem.projection)
        # Unmoved, the floods are mapped along the source lines (STRM_File_Clean), so they are burned there too
        if configs.relocate_misplaced_streams and configs.move_stream_network_to_thalweg:
            source_gdf = relocate_misplaced_streams(workspace, source_gdf, assigned_dem)
        else:
            workspace.relocated_StrmShp.unlink(missing_ok=True)
        _, dem_for_conflation = smooth_and_burn_dem(
            workspace, 
            source_gdf, 
            lakes,
            configs.stream_id_field, 
            configs.downstream_id_field
        )
        dem_for_conflation_path = workspace.fixed_dem
    elif should_move_streams:
        if configs.burn_streams:
            dem_for_conflation = Raster(workspace.fixed_dem).read_array()
            dem_for_conflation_path = workspace.fixed_dem
        else:
            dem_for_conflation = assigned_dem.read_array()
            dem_for_conflation_path = workspace.assigned_dem

        if configs.burn_streams and workspace.relocated_StrmShp.exists():
            # The DEM was burned on an earlier run, along these lines
            source_gdf = Vector(workspace.relocated_StrmShp, not configs.parallel).to_geopandas() \
                .to_crs(assigned_dem.projection)
        elif not configs.burn_streams and configs.relocate_misplaced_streams:
            # Unburned, the DEM's flow directions are the ones the streams are moved onto too
            derive_flow_directions(workspace, dem_for_conflation_path)
            flow_directions_ready = True
            source_gdf = Vector(workspace.DEM_StrmShp, not configs.parallel).to_geopandas().to_crs(assigned_dem.projection)
            source_gdf = relocate_misplaced_streams(workspace, source_gdf, assigned_dem, derived=True)
        else:
            source_gdf = Vector(workspace.DEM_StrmShp, not configs.parallel).to_geopandas().to_crs(assigned_dem.projection)

    if should_move_streams:
        if not flow_directions_ready:
            derive_flow_directions(workspace, dem_for_conflation_path)
        streams_gdf, final_streams = move_streams_onto_flow_directions(workspace, source_gdf, dem_for_conflation,
                                                                        assigned_dem)
        if streams_gdf.empty:
            if configs.raise_errors_if_nothing_in_domain:
                raise ValueError("No stream geometries remain after conflation.")
            else:
                workspace.new_StrmShp_matched.unlink(missing_ok=True)
                workspace.DEM_StrmShp = workspace.new_StrmShp_matched
                return None

        kwargs = {'index': False}
        if workspace.new_StrmShp_matched.suffix.lower().endswith('.parquet'):
            kwargs['compression'] = 'brotli'
            kwargs['write_covering_bbox'] = True
            kwargs['geometry_encoding'] = 'geoarrow'

        workspace.new_StrmShp_matched.parent.mkdir(parents=True, exist_ok=True)
        Vector.save_any_geom(streams_gdf, workspace.new_StrmShp_matched, **kwargs)
        if configs.minimize_output_files:
            workspace.DEM_StrmShp.unlink()
            workspace.relocated_StrmShp.unlink(missing_ok=True)
        _write_stream_raster(workspace.new_stream_raster, final_streams, assigned_dem)

        # The water mask is the moved streams and the water at the river's level beside them (see water_mask). ARC
        # carves bathymetry into all of it, so a line burned or drawn where no reach was moved onto (a relocated
        # line the water doesn't take, a reach left out) would be carved, and flood: it is no part of the mask.
        channel_mask = water_mask(workspace, assigned_dem, final_streams,
                                  reach_areas_km2(streams_gdf, configs)).astype(np.uint8)

        workspace.bathy_water_mask.parent.mkdir(parents=True, exist_ok=True)
        make_channel_mask(channel_mask, str(dem_for_conflation_path), str(workspace.bathy_water_mask), configs.compression)

    if configs.mapper.is_curve2flood_fldpln_mapper():
        if not should_move_streams:
            # Nothing was regenerated this run, so read back whichever stream raster the
            # rest of the pipeline is pointed at (see params["Stream_File"] in tasks/configs.py).
            if configs.move_stream_network_to_thalweg:
                final_streams = Raster(workspace.new_stream_raster).read_array()
            else:
                final_streams = Raster(workspace.STRM_File_Clean).read_array()

        _create_stream_info_table(
            workspace.stream_info_file,
            streams_array=final_streams,
            flow_direction_file=workspace.flowdir,
        )

def move_streams_onto_flow_directions(
        workspace: Workspace,
        source_gdf: gpd.GeoDataFrame,
        dem: np.ndarray,
        dem_raster: Raster) -> tuple[gpd.GeoDataFrame, np.ndarray]:
    """
    Move the source streams onto workspace.flowdir's D8 paths (see stream_conflation), and return them with
    their stream raster.

    Only flow paths draining at least new_strm_threshold_km2 are matched onto. For FLDPLN, which doesn't map
    lakes, lakes cut the network: no stream runs inside one. With the other mappers, reaches touching a lake are
    dropped unless they carry a river through it.
    """
    configs = workspace.configs
    flowdir_ds: gdal.Dataset = gdal.Open(str(workspace.flowdir))
    flowdir = flowdir_ds.ReadAsArray()
    flowdir_ds = None

    valid = np.isfinite(dem) & (dem != 0)
    if dem_raster.nodata_value is not None:
        valid &= dem != dem_raster.nodata_value
    lakes = load_lake_array(workspace, dem_raster)
    cut_lakes = lakes is not None and configs.mapper.is_curve2flood_fldpln_mapper()

    id_col, ds_col = configs.stream_id_field, configs.downstream_id_field
    source_gdf, area_col = _with_area_km2(source_gdf, configs)
    order_col = configs.StrmOrder_Field if configs.StrmOrder_Field in source_gdf.columns else None

    result = conflate_to_flow_directions(
        source_gdf, flowdir, dem_raster.geotransform, dem_raster.projection, id_col, ds_col,
        order_col=order_col, area_col=area_col, lake_mask=lakes if cut_lakes else None, valid_mask=valid,
        settings=ConflationSettings(min_channel_area_km2=configs.new_strm_threshold_km2),
    )
    streams, raster = result.streams, result.raster
    if area_col == "__area_km2":
        streams = streams.drop(columns=[area_col])

    if lakes is not None and not cut_lakes:
        stuck = _reaches_stuck_in_lakes(streams, raster, lakes, id_col, ds_col)
        streams, raster = _drop_reaches(streams, raster, stuck, id_col, ds_col)

    if configs.drop_multilinestrings:
        streams, raster = _drop_reaches(streams, raster, _multipart_reaches(workspace, source_gdf), id_col, ds_col)

    if configs.StrmOrder_Field and (configs.StrmOrder_Lower is not None or configs.StrmOrder_Upper is not None) and not configs.mapper.is_curve2flood_fldpln_mapper():
        kept = _filter_streams_by_stream_order(streams.copy(), configs.StrmOrder_Field, configs.StrmOrder_Lower, configs.StrmOrder_Upper)
        streams, raster = _drop_reaches(streams, raster, set(streams[id_col]) - set(kept[id_col]), id_col, ds_col)

    return streams, raster

def _multipart_reaches(workspace: Workspace, source_gdf: gpd.GeoDataFrame) -> set:
    """
    The ids of the reaches whose line is in pieces that can't be joined into one, in the streams as they were before
    relocation (workspace.DEM_StrmShp, or else ``source_gdf``): relocation redraws such a reach as one line, the path
    of its main stem. Decided after a line_merge, not from the stored geometry type: a GeoPackage stores every
    LineString as a MultiLineString.
    """
    if workspace.DEM_StrmShp.exists():
        source_gdf = Vector(workspace.DEM_StrmShp, not workspace.configs.parallel).to_geopandas()
    merged = source_gdf.geometry.line_merge()
    return set(source_gdf.loc[(merged.geom_type == 'MultiLineString').to_numpy(), workspace.configs.stream_id_field])


def _with_area_km2(source_gdf: gpd.GeoDataFrame, configs) -> tuple[gpd.GeoDataFrame, str | None]:
    """
    ``source_gdf`` and the name of its drainage area column in km2: area_km2_field, or area_m2_field converted
    into a '__area_km2' column. None when it has neither.
    """
    if configs.area_km2_field and configs.area_km2_field in source_gdf.columns:
        return source_gdf, configs.area_km2_field
    if configs.area_m2_field and configs.area_m2_field in source_gdf.columns:
        area = pd.to_numeric(source_gdf[configs.area_m2_field], errors="coerce") / 1e6
        return source_gdf.assign(__area_km2=area), "__area_km2"
    return source_gdf, None

def relocate_misplaced_streams(workspace: Workspace, source_gdf: gpd.GeoDataFrame, dem_raster: Raster,
                               derived: bool = False) -> gpd.GeoDataFrame:
    """
    Move the reaches the unburned DEM shows are off their river onto its D8 flow paths (see stream_relocation), and
    save the network to workspace.relocated_StrmShp when any line changed. Where the DEM dams a river, the source
    lines are burned through first (see breach_dammed_rivers).

    When ``derived``, workspace.filled_dem and workspace.flowdir are the unburned DEM's, and are used unless it
    dams a river; otherwise the flow directions are derived here.
    """
    configs = workspace.configs
    workspace.relocated_StrmShp.unlink(missing_ok=True)
    gdf, area_col = _with_area_km2(source_gdf, configs)
    if area_col is None:
        LOG.warning("Misplaced streams weren't relocated: the streams have no drainage area "
                    "(area_km2_field or area_m2_field).")
        return source_gdf

    dem = dem_raster.read_array()
    valid = np.isfinite(dem) & (dem != 0)
    if dem_raster.nodata_value is not None:
        valid &= dem != dem_raster.nodata_value
    breached, _ = breach_dammed_rivers(gdf, dem, dem_raster.geotransform, dem_raster.projection,
                                       configs.stream_id_field, configs.downstream_id_field, valid_mask=valid)
    ground = dem if breached is None else breached
    if derived and breached is None:
        flowdir_array, filled_array = Raster(workspace.flowdir).read_array(), Raster(workspace.filled_dem).read_array()
    else:
        filled_array, flowdir_array = fill_and_flow_directions(ground, valid, dem_raster.geotransform,
                                                               dem_raster.projection)
    filled_array = np.where(valid, filled_array, np.nan)
    streams = relocate_misplaced_reaches(gdf, flowdir_array, filled_array, ground, dem_raster.geotransform,
                                         dem_raster.projection, configs.stream_id_field, configs.downstream_id_field,
                                         area_col, valid_mask=valid)
    streams = streams.drop(columns=["__area_km2"], errors="ignore")
    changed = ~streams.geometry.geom_equals(source_gdf.geometry).to_numpy()
    if not changed.any():
        return source_gdf
    kwargs = {'index': False}
    if workspace.relocated_StrmShp.suffix.lower() == '.parquet':
        kwargs.update(compression='brotli', write_covering_bbox=True, geometry_encoding='geoarrow')
    workspace.relocated_StrmShp.parent.mkdir(parents=True, exist_ok=True)
    Vector.save_any_geom(streams, workspace.relocated_StrmShp, **kwargs)

    _replace_lines_in_water_mask(workspace, dem_raster, streams[changed], streams.loc[changed, configs.stream_id_field])
    return streams

def _replace_lines_in_water_mask(workspace: Workspace, dem_raster: Raster, redrawn: gpd.GeoDataFrame,
                                 ids: pd.Series) -> None:
    """
    Take the ``redrawn`` reaches' old lines out of the water mask, except where the land cover has water, and put
    their new lines in: the burn smooths the channel along the mask's cells, which should be the lines it burns.
    """
    if not workspace.bathy_water_mask.exists() or not workspace.STRM_File_Clean.exists():
        return
    mask = Raster(workspace.bathy_water_mask).read_array().astype(bool)
    old = np.isin(Raster(workspace.STRM_File_Clean).read_array(), pd.to_numeric(ids).to_numpy())
    if workspace.LAND_File.exists():
        old &= Raster(workspace.LAND_File).read_array() != workspace.configs.land_watervalue

    new_ds: gdal.Dataset = gdal.GetDriverByName("MEM").Create("", dem_raster.ncols, dem_raster.nrows, 1, gdal.GDT_Byte)
    new_ds.SetGeoTransform(dem_raster.geotransform)
    new_ds.SetProjection(dem_raster.projection)
    ogr_ds, layer = _in_memory_ogr_layer(redrawn, dem_raster.projection)
    gdal.RasterizeLayer(new_ds, [1], layer, burn_values=[1])
    new = new_ds.ReadAsArray().astype(bool)
    ogr_ds = new_ds = None

    mask = (mask & ~old) | new
    make_channel_mask(mask.astype(np.uint8), str(workspace.assigned_dem), str(workspace.bathy_water_mask),
                      workspace.configs.compression)

def _reaches_stuck_in_lakes(streams: gpd.GeoDataFrame, raster: np.ndarray, lakes: np.ndarray, id_col: str,
                            ds_col: str) -> set:
    """
    The reaches that touch a lake, except those that carry a river through it: ones with a reach somewhere above
    them and a reach somewhere below them that touch no lake.
    """
    in_lake = set(np.unique(raster[lakes & (raster > 0)]).tolist())
    if not in_lake:
        return set()
    graph = nx.DiGraph()
    graph.add_nodes_from(streams[id_col])
    graph.add_edges_from((r, d) for r, d in zip(streams[id_col], streams[ds_col]) if d in graph and d != r)
    stuck = set()
    for reach in in_lake:
        if not (nx.descendants(graph, reach) - in_lake) or not (nx.ancestors(graph, reach) - in_lake):
            stuck.add(reach)
    return stuck

def _drop_reaches(streams: gpd.GeoDataFrame, raster: np.ndarray, reaches: set, id_col: str,
                  ds_col: str) -> tuple[gpd.GeoDataFrame, np.ndarray]:
    """Remove ``reaches`` from the conflated network and its raster; reaches that drained into them become outlets."""
    if not reaches:
        return streams, raster
    raster = np.where(np.isin(raster, list(reaches)), 0, raster).astype(raster.dtype, copy=False)
    streams = streams[~streams[id_col].isin(reaches)].copy()
    streams.loc[streams[ds_col].isin(reaches), ds_col] = -1
    return streams, raster

def _write_stream_raster(path: Path, streams: np.ndarray, reference: Raster) -> None:
    """Write the conflated stream ids as an Int32 GeoTIFF on ``reference``'s grid, 0 off the streams."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    driver: gdal.Driver = gdal.GetDriverByName("GTiff")
    if path.exists():
        driver.Delete(str(path))
    stream_ds: gdal.Dataset = driver.Create(str(path), streams.shape[1], streams.shape[0], 1, gdal.GDT_Int32,
                                            ['COMPRESS=DEFLATE', 'PREDICTOR=2'])
    stream_ds.SetGeoTransform(reference.geotransform)
    stream_ds.SetProjection(reference.projection)
    stream_ds.GetRasterBand(1).WriteArray(streams.astype(np.int32, copy=False))
    # Closed explicitly: a numba compile earlier in the run can keep this frame, and so the dataset, alive until the
    # next garbage collection, and until then the raster isn't written
    stream_ds.Close()

def smooth_and_burn_dem(
        workspace: Workspace, 
        source_gdf: gpd.GeoDataFrame,
        lakes: np.ndarray | None = None, 
        id_col: str = 'LINKNO',
        ds_col: str = 'DSLINKNO') -> tuple[np.ndarray, np.ndarray]:
    """
    Burn ``source_gdf``'s streams into the DEM, smooth the burned channels, and write workspace.fixed_dem.
    Returns the channel mask and the burned DEM. With use_dem_derived_channel_mask, the water the DEM flattens along
    the streams (see channel_detection) is smoothed as water, and the burn leaves the lines on it (see
    burn_streams_into_dem).
    """
    dem_ds: gdal.Dataset = gdal.Open(workspace.assigned_dem)
    dem = dem_ds.ReadAsArray()
    nodata_value = dem_ds.GetRasterBand(1).GetNoDataValue()
    if nodata_value is None:
        nodata_value = -9999
    nan_mask = np.isnan(dem)
    dem[nan_mask] = nodata_value

    ocean_mask = (dem == 0)
    dem[ocean_mask] = nodata_value

    configs = workspace.configs
    area_col = next((c for c in (configs.area_km2_field, configs.area_m2_field) if c and c in source_gdf.columns), None)
    streams = _get_stream_raster(dem_ds, source_gdf, id_col, configs.StrmOrder_Field, area_col)

    channel_mask = Raster(workspace.bathy_water_mask).read_array()
    leveled = None
    if workspace.configs.use_dem_derived_channel_mask:
        # The water the DEM flattens along the streams (see channel_detection). Not the low ground by them as well:
        # measured from the lines being burned, that would take in every line the burn has to lower
        leveled = dem_channel(dem_ds, dem, (dem != nodata_value) & np.isfinite(dem), streams,
                              reach_areas_km2(source_gdf, configs)).flat_water
        channel_mask = channel_mask.astype(bool) | leveled

    original = dem.copy()
    paths = []
    kept = np.zeros(dem.shape, bool)
    channel_mask = burn_streams_into_dem(channel_mask, dem, streams, dem_ds, source_gdf, lakes, id_col, ds_col,
                                         nodata_value, order_col=configs.StrmOrder_Field, area_col=area_col,
                                         paths=paths, leveled=leveled, kept=kept)
    del leveled
    dem = smooth_burned_dem(dem, channel_mask, streams, original, nodata_value, paths=paths, fixed=kept)

    dem[dem < -1000] = nodata_value # Remove any DEM values that are less than -1000 m, since these are likely to be erroneous and will cause problems with the floodplain mapping.
    output_ds = gdal.GetDriverByName('GTiff').Create(workspace.fixed_dem, dem_ds.RasterXSize, dem_ds.RasterYSize, 1, gdal.GDT_Float32, options=[f'COMPRESS={workspace.configs.compression}'])
    output_ds.WriteArray(dem)
    output_ds.SetGeoTransform(dem_ds.GetGeoTransform())
    output_ds.SetProjection(dem_ds.GetProjection())
    output_ds.GetRasterBand(1).SetNoDataValue(nodata_value)
    output_ds = None

    return channel_mask, dem

@njit(cache=True)
def _label_terraces(values: np.ndarray, use: np.ndarray, n_use: int) -> tuple[np.ndarray, np.ndarray, int]:
    """
    Number the 8-connected groups of ``use`` cells that share one value (0 elsewhere), and list their cells as flat
    indices, group by group.
    """
    nrows, ncols = values.shape
    labels = np.zeros((nrows, ncols), np.int32)
    cells = np.empty(max(n_use, 1), np.int64)
    n = 0
    count = 0
    for r0 in range(nrows):
        for c0 in range(ncols):
            if not use[r0, c0] or labels[r0, c0]:
                continue
            count += 1
            z = values[r0, c0]
            labels[r0, c0] = count
            cells[n] = r0 * ncols + c0
            head = n
            n += 1
            while head < n:
                r, c = cells[head] // ncols, cells[head] % ncols
                head += 1
                for nr in range(max(r - 1, 0), min(r + 2, nrows)):
                    for nc in range(max(c - 1, 0), min(c + 2, ncols)):
                        if use[nr, nc] and labels[nr, nc] == 0 and values[nr, nc] == z:
                            labels[nr, nc] = count
                            cells[n] = nr * ncols + nc
                            n += 1
    return labels, cells[:n], count

@njit(cache=True)
def _distance_to_step(values: np.ndarray, labels: np.ndarray, cells: np.ndarray, index: np.ndarray, max_step: float,
                      up: bool, along: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    For each of ``cells`` (``index`` gives each cell's place among them, -1 elsewhere): how many cells it is, within
    its terrace, from the nearest of its terrace's cells next to a terrace higher (``up``) or lower than it by at most
    ``max_step``, and that terrace's elevation. -1 where there is none. A cell's ``along`` (one per cell, NaN for
    none), the elevation of the cell above it (``up``) or below it along a path, is such a terrace whatever the step,
    and is taken over the cells around it.
    """
    nrows, ncols = values.shape
    n = cells.size
    dist = np.full(n, -1, np.int32)
    other = np.zeros(n, np.float32)
    queue = np.empty(max(n, 1), np.int64)
    tail = 0
    for i in range(n):
        r, c = cells[i] // ncols, cells[i] % ncols
        label, z = labels[r, c], values[r, c]
        best = np.inf
        step = along[i] - z if up else z - along[i]
        if step > 0:  # False for NaN
            best = step
            other[i] = along[i]
        else:
            for nr in range(max(r - 1, 0), min(r + 2, nrows)):
                for nc in range(max(c - 1, 0), min(c + 2, ncols)):
                    if labels[nr, nc] == 0 or labels[nr, nc] == label:
                        continue
                    step = values[nr, nc] - z if up else z - values[nr, nc]
                    if 0 < step <= max_step and step < best:
                        best = step
                        other[i] = values[nr, nc]
        if best < np.inf:
            dist[i] = 0
            queue[tail] = i
            tail += 1
    head = 0
    while head < tail:
        i = queue[head]
        head += 1
        r, c = cells[i] // ncols, cells[i] % ncols
        for nr in range(max(r - 1, 0), min(r + 2, nrows)):
            for nc in range(max(c - 1, 0), min(c + 2, ncols)):
                k = index[nr, nc]
                if k >= 0 and dist[k] < 0 and labels[nr, nc] == labels[r, c]:
                    dist[k] = dist[i] + 1
                    other[k] = other[i]
                    queue[tail] = k
                    tail += 1
    return dist, other

@njit(cache=True)
def _distance_to_bank(cells: np.ndarray, index: np.ndarray) -> np.ndarray:
    """For each of ``cells``: how many cells it is from the nearest cell that isn't one of them (1 next to one)."""
    nrows, ncols = index.shape
    n = cells.size
    dist = np.zeros(n, np.int32)
    queue = np.empty(max(n, 1), np.int64)
    tail = 0
    for i in range(n):
        r, c = cells[i] // ncols, cells[i] % ncols
        for nr in range(max(r - 1, 0), min(r + 2, nrows)):
            for nc in range(max(c - 1, 0), min(c + 2, ncols)):
                if index[nr, nc] < 0 and dist[i] == 0:
                    dist[i] = 1
                    queue[tail] = i
                    tail += 1
    head = 0
    while head < tail:
        i = queue[head]
        head += 1
        r, c = cells[i] // ncols, cells[i] % ncols
        for nr in range(max(r - 1, 0), min(r + 2, nrows)):
            for nc in range(max(c - 1, 0), min(c + 2, ncols)):
                k = index[nr, nc]
                if k >= 0 and dist[k] == 0:
                    dist[k] = dist[i] + 1
                    queue[tail] = k
                    tail += 1
    return dist

def _ramp_terraces(values: np.ndarray, use: np.ndarray, max_step: float, min_cells: int = 1,
                   above: tuple[np.ndarray, np.ndarray] | None = None,
                   below: tuple[np.ndarray, np.ndarray] | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Turn the terraces among ``use`` cells, groups of cells sharing one elevation, into ramps. Where a terrace meets
    one up to ``max_step`` higher or lower, it takes the mean of the two; in its middle, farthest from both, it keeps
    its own elevation; in between, it is interpolated by distance within the terrace. A terrace meeting only higher
    or only lower ones keeps its elevation where it is farthest from them.

    ``above`` and ``below``, if given, are cells (flat indices) and the elevation of the cell above or below each
    along a path (see _path_steps). A terrace meets the next one there whatever the step, rather than among the
    cells around it; with a ``max_step`` of 0, only there.

    Returns the cells of the terraces of at least ``min_cells`` cells (flat indices), their new elevations, and how
    near each is to the middle of its terrace, away from the cells off the terraces: 0 next to one, 1 farthest.
    """
    labels, cells, count = _label_terraces(values, use, int(np.count_nonzero(use)))
    terrace = labels.ravel()[cells]
    if min_cells > 1:
        small = np.bincount(terrace, minlength=count + 1)[terrace] < min_cells
        labels.ravel()[cells[small]] = 0
        cells, terrace = cells[~small], terrace[~small]
    if cells.size == 0:
        return cells, np.empty(0), np.empty(0)
    index = np.full(values.shape, -1, np.int32)
    index.ravel()[cells] = np.arange(cells.size, dtype=np.int32)
    along = []
    for steps, nearest in ((above, np.fmin), (below, np.fmax)):
        z = np.full(cells.size, np.nan)
        if steps is not None:
            k = index.ravel()[steps[0]]
            nearest.at(z, k[k >= 0], steps[1][k >= 0])  # the smallest step, where a cell has several
        along.append(z)
    up, higher = _distance_to_step(values, labels, cells, index, max_step, True, along[0])
    down, lower = _distance_to_step(values, labels, cells, index, max_step, False, along[1])
    bank = _distance_to_bank(cells, index).astype(np.float64)
    del labels, index

    z = values.ravel()[cells].astype(np.float64)
    higher, lower = higher.astype(np.float64), lower.astype(np.float64)
    has_up, has_down = up >= 0, down >= 0
    # Measured from the edge between two terraces, so the cells on either side of it are equally far from it
    up, down = up + 0.5, down + 0.5
    farthest_up = np.zeros(count + 1)
    np.maximum.at(farthest_up, terrace[has_up], up[has_up])
    farthest_down = np.zeros(count + 1)
    np.maximum.at(farthest_down, terrace[has_down], down[has_down])

    new = z.copy()
    k = has_up & has_down & (up <= down)
    new[k] += (higher[k] - z[k]) / 2 * (down[k] - up[k]) / (down[k] + up[k])
    k = has_up & has_down & (up > down)
    new[k] += (lower[k] - z[k]) / 2 * (up[k] - down[k]) / (down[k] + up[k])
    k = has_up & ~has_down
    new[k] += (higher[k] - z[k]) / 2 * (1 - up[k] / farthest_up[terrace[k]])
    k = has_down & ~has_up
    new[k] += (lower[k] - z[k]) / 2 * (1 - down[k] / farthest_down[terrace[k]])

    widest = np.zeros(count + 1)
    np.maximum.at(widest, terrace, bank)
    span = widest[terrace] - 1
    centre = np.divide(bank - 1, span, out=np.zeros_like(bank), where=span > 0)
    return cells, new, centre

def _path_steps(values: np.ndarray, paths: list) -> tuple[tuple[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]:
    """
    Where ``values`` (the burned DEM) falls from one cell to the next down ``paths`` (as smooth_burned_dem takes
    them, a path's start cells above its first), for _ramp_terraces: the cells below each drop with the elevation
    above them (``above``), and the cells above each drop with the elevation below them (``below``).

    The cells below a drop are only those on their own path, the last one burned through them: a tributary's line
    ends on the river's, and the river shouldn't ramp up toward it there.
    """
    flat = values.ravel()
    own = np.full(flat.size, -1, np.int32)
    for i, (_, cells) in enumerate(paths):
        own[cells] = i
    tops, under_tops, bottoms, over_bottoms = [], [], [], []
    for i, (starts, cells) in enumerate(paths):
        upper = np.concatenate([starts, cells[:-1]])
        lower = np.concatenate([np.full(starts.size, cells[0]), cells[1:]])
        falls = flat[upper] > flat[lower]
        upper, lower = upper[falls], lower[falls]
        tops.append(upper)
        under_tops.append(flat[lower])
        mine = own[lower] == i
        bottoms.append(lower[mine])
        over_bottoms.append(flat[upper[mine]])
    above = (np.concatenate(bottoms), np.concatenate(over_bottoms).astype(np.float64))
    below = (np.concatenate(tops), np.concatenate(under_tops).astype(np.float64))
    return above, below

def smooth_burned_dem(
        dem: np.ndarray,
        mask: np.ndarray | None = None,
        streams: np.ndarray | None = None,
        original: np.ndarray | None = None,
        nodata_value: float = -9999,
        step: float = 0.5,
        centre_drop: float = 0.1,
        min_water_cells: int = 6,
        paths: list | None = None,
        fixed: np.ndarray | None = None) -> np.ndarray:
    """
    Smooth the stair steps out of the channels of a burned DEM, and return it.

    Some DEMs, like FABDEM, flatten water into terraces at whole multiples of ``step``, and the burn lowers the stream
    lines in such steps too. Stair steps throw off the stream slopes ARC measures, so each terrace becomes a ramp:
    where it meets the next it takes their mean, and in its middle it keeps its own elevation (see _ramp_terraces).

    The water is the ``mask``'s cells at whole multiples of ``step`` in ``original``, the DEM before the burn, so the
    cells the burn lowered or raised in it are smoothed with the water around them. Its middle, farthest from the
    banks, is lowered by up to ``centre_drop`` for the flow directions. The burned lines (``streams``) are smoothed
    along themselves, and in the water kept as low as its middle, or lower where the burn cut deeper. A bank no
    higher than a smoothed cell beside it is raised 0.5 m above it.

    Given ``paths``, a line's terraces meet where it drops down them, by however much. The burn keeps a line no
    higher than the lowest cell beside it upstream, so a line can run level for a whole reach and then drop metres
    where the ground does; ramped only toward steps up to ``step``, as the water is, such a terrace stays level.

    The water's ramps can rise along a line the burn kept level, leaving pools along it, so the lines are then made
    never to rise along ``paths``, the lines burn_streams_into_dem burned, in the order it burned them, and to keep
    falling as their own ramps fall. Kept level instead, a line across a pool FABDEM flattened can drain either way
    along itself, and the flow directions take a river back out over the DEM's edge it came in over. The cells the
    burn left (``fixed``, see burn_streams_into_dem) are smoothed as the water around them, and aren't lowered there
    either: a line goes on from where they are.
    """
    dem = dem.astype(np.float32, copy=False)
    original = dem if original is None else original.astype(np.float32, copy=False)
    if mask is None:
        mask = (original % step) == 0
    mask = mask.astype(bool)
    valid = np.isfinite(dem) & (dem != nodata_value)

    water_cells, water, centre = _ramp_terraces(original, mask & valid & (original % step == 0), step, min_water_cells)
    if streams is None:
        line_cells, line = np.empty(0, np.int64), np.empty(0)
    else:
        # A line the burn left on the water (``fixed``) is smoothed as the water around it, not lowered to its middle
        lines = (streams > 0) & valid & (dem % step == 0)
        if fixed is not None:
            lines &= ~np.asarray(fixed, bool)
        if paths:
            # Only the steps down the lines' own paths: lines that touch elsewhere (a tributary beside the river, a
            # bend come back on itself) aren't the next step down either one
            above, below = _path_steps(dem, paths)
            line_cells, line, _ = _ramp_terraces(dem, lines, 0, above=above, below=below)
        else:
            line_cells, line, _ = _ramp_terraces(dem, lines, step)

    out = dem.copy()
    out.ravel()[water_cells] = water - centre_drop * centre
    in_water = np.full(dem.size, -1, np.int32)
    in_water[water_cells] = np.arange(water_cells.size, dtype=np.int32)
    k = in_water[line_cells]
    del in_water
    wet = k >= 0
    ramp = line.copy()
    line[wet] = np.minimum(line[wet], water[k[wet]] - centre_drop)
    out.ravel()[line_cells] = line

    changed = np.zeros(dem.shape, bool)
    changed.ravel()[water_cells] = True
    changed.ravel()[line_cells] = True
    beside = maximum_filter(np.where(changed, out, -np.inf).astype(np.float32), size=3)
    banks = ~mask & ~changed & valid & (out <= beside)
    out[banks] = beside[banks] + 0.5
    if fixed is not None:
        # Where the burn lowered one of them all the same (a line's end at the DEM's edge), it stays as low
        lowered = np.asarray(fixed, bool) & (dem < original)
        out[lowered] = np.minimum(out[lowered], dem[lowered])

    if paths:
        ramps = out.copy()
        ramps.ravel()[line_cells] = ramp
        starts = [start for start, _ in paths]
        cells = [path for _, path in paths]
        _never_rise(out.ravel(), ramps.ravel(), np.concatenate(starts), np.cumsum([s.size for s in starts]),
                    np.concatenate(cells), np.cumsum([c.size for c in cells]),
                    np.zeros(0, np.bool_) if fixed is None else np.ascontiguousarray(fixed, np.bool_).ravel())
    return out

@njit(cache=True)
def _never_rise(values: np.ndarray, ramps: np.ndarray, starts: np.ndarray, start_ends: np.ndarray,
                cells: np.ndarray, cell_ends: np.ndarray, fixed: np.ndarray) -> None:
    """
    Lower ``values`` (flat) along each path in turn, so that none rises, and each keeps falling where its ramp
    (``ramps``, at least ``values``) falls: path i runs down cells[cell_ends[i - 1]:cell_ends[i]], starting no higher
    than its start cells, starts[start_ends[i - 1]:start_ends[i]]. A path stays as far below its ramp as it has had
    to be anywhere above, its start included. The ``fixed`` cells (flat; none if empty) aren't lowered, and the path
    goes on from them as from its start.
    """
    first_start, first_cell = 0, 0
    for i in range(cell_ends.size):
        level = np.inf
        for k in range(first_start, start_ends[i]):
            level = min(level, values[starts[k]])
        below = ramps[cells[first_cell]] - min(values[cells[first_cell]], level)
        for k in range(first_cell, cell_ends[i]):
            cell = cells[k]
            if fixed.size and fixed[cell]:
                level = values[cell]
                below = ramps[cell] - values[cell]
                continue
            below = max(below, ramps[cell] - values[cell])
            level = min(level, ramps[cell] - below)
            values[cell] = level
        first_start, first_cell = start_ends[i], cell_ends[i]

def _get_stream_raster(dem_ds: gdal.Dataset, streams_gdf: gpd.GeoDataFrame, id_col: str = 'LINKNO',
                       order_col: str | None = None, area_col: str | None = None) -> np.ndarray:
    """
    Rasterize the stream ids onto the DEM's grid.

    A cell where reaches meet holds the last one drawn, and is burned only as deep as that reach. The larger
    reaches (by stream order, then drainage area) are drawn last, so a river keeps the cells it shares with its
    tributaries: a tributary's cell would dam the river, and filling would turn the channel above into a pool.
    """
    keys = [c for c in (order_col, area_col) if c and c in streams_gdf.columns]
    if keys:
        streams_gdf = streams_gdf.sort_values(keys, kind="stable", na_position="first")
    mem_ds = gdal.GetDriverByName('MEM').Create('', dem_ds.RasterXSize, dem_ds.RasterYSize, 1, gdal.GDT_Int32)
    mem_ds.SetGeoTransform(dem_ds.GetGeoTransform())
    mem_ds.SetProjection(dem_ds.GetProjection())
    ogr_ds, stream_layer = _in_memory_ogr_layer(streams_gdf, dem_ds.GetProjection(), values=streams_gdf[id_col])
    gdal.RasterizeLayer(mem_ds, [1], stream_layer, options=["ATTRIBUTE=value"])
    ogr_ds = None
    mem_ds.FlushCache()
    streams = mem_ds.ReadAsArray()

    return streams

def burn_streams_into_dem(
        channel_mask: np.ndarray,
        dem: np.ndarray,
        streams: np.ndarray,
        dem_ds: gdal.Dataset,
        streams_gdf: gpd.GeoDataFrame,
        lakes: np.ndarray | None,
        id_col: str = 'LINKNO',
        ds_col: str = 'DSLINKNO',
        nodata_value: float = -9999,
        min_feature_size: int = 5,
        order_col: str | None = None,
        area_col: str | None = None,
        paths: list | None = None,
        leveled: np.ndarray | None = None,
        kept: np.ndarray | None = None,
):
    """
    Burn every reach into ``dem`` (in place), upstream reaches first, in 0.5 m steps that never rise downstream,
    each starting no higher than the largest of the reaches flowing into it, by ``order_col`` and then ``area_col``
    (all of them, without either). Returns the channel mask without its specks (``min_feature_size`` cells or
    fewer), lakes, ocean and nodata. Appends each burned line to ``paths``, if given (see burn_linestring).

    A line isn't lowered on ``leveled``, if given, the water the DEM flattens (see channel_detection): the line goes
    on down from such a cell at its level. Except at a line's end on the DEM's edge, which is lowered so the flow leaves there. The cells left so are
    marked in ``kept``, if given.
    """
    G = nx.from_pandas_edgelist(
        streams_gdf[streams_gdf[ds_col] > 0],
        source=id_col,
        target=ds_col,
        create_using=nx.DiGraph
    )
    # A reach with nothing flowing in or out of it is burned too
    G.add_nodes_from(streams_gdf[id_col])
    # Remove any cells that have less than 5 connected neighbors in the mask, as they are likely to be noise or small artifacts
    structure = np.ones((3, 3), dtype=int)
    labels, _ = label(channel_mask, structure=structure)

    counts = np.bincount(labels.ravel())

    keep = counts > min_feature_size
    keep[0] = False
    channel_mask = keep[labels]

    # Mask out lakes
    if lakes is not None:
        channel_mask = channel_mask.astype(bool) & ~lakes

    # Mask out ocean (where elevation == 0)
    ocean_mask = (dem == 0)
    channel_mask &= ~ocean_mask

    # Mask out dem's no data value
    if nodata_value is not None:
        channel_mask &= (dem != nodata_value)

    # Buffer the stream raster as a new mask
    channel_border = binary_dilation(streams > 0, structure=structure).astype(np.uint8, copy=False)
    channel_border &= (streams == 0) & (dem > -9998)

    # Burn upstream reaches first. A reach starts no higher than the reaches flowing into it (see
    # _largest_upstream_reaches), which only holds if they were burned before it; in table order a river could be
    # burned before its tributaries had been, leaving a step up where they meet.
    try:
        burn_order = {stream_id: i for i, stream_id in enumerate(nx.topological_sort(G))}
        streams_gdf = streams_gdf.iloc[np.argsort(streams_gdf[id_col].map(burn_order).to_numpy(), kind="stable")]
    except nx.NetworkXUnfeasible:
        LOG.warning("The source stream network has a cycle, so its reaches are burned in table order.")
    streams_gdf = streams_gdf.set_index(id_col)

    masked_dem = np.where(channel_border, dem, np.inf)
    local_min = minimum_filter(masked_dem, size=3, mode="nearest")
    upstream = _largest_upstream_reaches(G, streams_gdf, order_col, area_col)

    # Traverse each segment. Identify which end is up and downstream.
    # Then, take the upstream value. If it is not a multiple of 0.5, lower it to so.
    # Go downstream. Each cell is a multiple of 0.5, no higher than upstream elevation or the current cell elevation.
    for row in streams_gdf.itertuples():
        geom = row.geometry
        stream_id = row.Index
        if geom.geom_type == "MultiLineString":
            lines = geom.geoms
        elif geom.geom_type == "LineString":
            lines = [geom]
        else:
            raise ValueError(f"Unsupported geometry type: {geom.geom_type}")
        for line in lines:
            path = burn_linestring(dem, streams, dem_ds, line, stream_id, G, streams_gdf, local_min, nodata_value,
                                   upstream.get(stream_id, []), leveled, kept)
            if path is not None and paths is not None:
                paths.append(path)

    return channel_mask

def _largest_upstream_reaches(G: nx.DiGraph, streams_gdf: gpd.GeoDataFrame, order_col: str | None = None,
                              area_col: str | None = None) -> dict:
    """
    Each reach's largest upstream reaches (``streams_gdf`` is indexed by stream id): those flowing into it of the
    highest ``order_col`` and then ``area_col``, as _get_stream_raster ranks them, or all of them without either.

    A reach starts no higher than these. Taking every reach flowing in let a tributary whose mouth lies below the
    river's water, as FABDEM's flattened rivers often have, lower the river below it for kilometres. Such a
    tributary is left in a pool at its mouth instead, as the river's backwater would leave it.
    """
    keys = [c for c in (order_col, area_col) if c and c in streams_gdf.columns]
    sizes = {}
    if keys:
        table = streams_gdf[keys].apply(pd.to_numeric, errors="coerce").fillna(-np.inf)
        table = table.groupby(level=0).max()
        sizes = dict(zip(table.index, table.itertuples(index=False, name=None)))
    upstream = {}
    for reach in G:
        reaches = list(G.predecessors(reach))
        if sizes and len(reaches) > 1:
            largest = max(sizes.get(r, ()) for r in reaches)
            reaches = [r for r in reaches if sizes.get(r, ()) == largest]
        upstream[reach] = reaches
    return upstream

@njit(cache=True, nogil=True)
def nearest_stream_raster_pixel(
    streams_array: np.ndarray,
    start_row: int,
    start_col: int,
    linkno: int
):
    """The cell of reach ``linkno`` nearest the start cell, within one cell of it, or the start cell if none."""
    arr = streams_array
    nrows, ncols = arr.shape

    best_dist = np.inf
    best_r = start_row
    best_c = start_col

    # manually unrolled 3x3 neighborhood
    for r in (start_row - 1, start_row, start_row + 1):
        if r < 0 or r >= nrows:
            continue
        for c in (start_col - 1, start_col, start_col + 1):
            if c < 0 or c >= ncols:
                continue

            if arr[r, c] == linkno:
                dr = c - start_col
                dc = r - start_row
                dist = dr * dr + dc * dc

                if dist < best_dist:
                    best_dist = dist
                    best_r = r
                    best_c = c

    return best_r, best_c

def burn_linestring(
        dem: np.ndarray, 
        streams: np.ndarray, 
        dem_ds: gdal.Dataset, 
        linestring: LineString, 
        linkno: int, 
        G: nx.DiGraph,
        streams_gdf: gpd.GeoDataFrame,
        local_min: np.ndarray,
        nodata_value: float,
        upstream_ids: list | None = None,
        leveled: np.ndarray | None = None,
        kept: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Burn one line of reach ``linkno`` into ``dem``, from its upstream end down, starting no higher than the
    reaches ``upstream_ids`` (default: all those flowing into it). The upstream end is found from the stream
    network, or else from the elevations. Its cells on the ``leveled`` water are left, and marked in ``kept``.

    Returns the cells of ``upstream_ids`` its start was kept below and the cells it burned, from upstream down,
    as flat indices, or None if it burned nothing.
    """
    if linkno not in G:
        return None  # Skip if the linkno is not in the graph
    
    inverse_transform = gdal.InvGeoTransform(dem_ds.GetGeoTransform())
    coords = np.asarray(linestring.coords)

    # Find the first and last points of the linestring in pixel coordinates
    x1, y1 = coords[0]
    col1, row1 = gdal.ApplyGeoTransform(inverse_transform, x1, y1)
    x2, y2 = coords[-1]
    col2, row2 = gdal.ApplyGeoTransform(inverse_transform, x2, y2)

    # Determine which end is upstream and which is downstream based on topolgy
    if not (0 <= row1 < dem.shape[0] and 0 <= col1 < dem.shape[1] and dem[math.floor(row1), math.floor(col1)] != nodata_value):
        # Traverse the linestring to find the first point that is within the DEM bounds
        for x, y in coords:
            col, row = gdal.ApplyGeoTransform(inverse_transform, x, y)
            if 0 <= row < dem.shape[0] and 0 <= col < dem.shape[1] and dem[math.floor(row), math.floor(col)] != nodata_value:
                col1, row1 = col, row
                break
        else:
            return None  # No valid point found within DEM bounds

    if not (0 <= row2 < dem.shape[0] and 0 <= col2 < dem.shape[1] and dem[math.floor(row2), math.floor(col2)] != nodata_value):
        # Traverse the linestring in reverse to find the last point that is within the DEM bounds
        for x, y in np.flipud(coords):
            col, row = gdal.ApplyGeoTransform(inverse_transform, x, y)
            if 0 <= row < dem.shape[0] and 0 <= col < dem.shape[1] and dem[math.floor(row), math.floor(col)] != nodata_value:
                col2, row2 = col, row
                break
        else:
            return None  # No valid point found within DEM bounds

    row1, col1 = math.floor(row1), math.floor(col1)
    row2, col2 = math.floor(row2), math.floor(col2)
    row1, col1 = nearest_stream_raster_pixel(streams, row1, col1, linkno)
    row2, col2 = nearest_stream_raster_pixel(streams, row2, col2, linkno)

    # Can't rely on elevation to determine upstream/downstream because the DEM may have been modified by previous
    # burns, and a flattened river's two ends are often level. Instead, use the topology of the stream network: the
    # upstream end is the one nearer the reaches flowing in, or else the one farther from the reach this one flows
    # into. Measured by distance rather than by whether an end exactly touches its neighbour, which a line that
    # stops a hair short of its junction never does.
    first_point, last_point = Point(x1, y1), Point(x2, y2)
    flip = None
    neighbours = [n for n in G.predecessors(linkno) if n in streams_gdf.index]
    if neighbours:
        geoms = [_first_geometry(streams_gdf, n) for n in neighbours]
        d_first = min(first_point.distance(g) for g in geoms)
        d_last = min(last_point.distance(g) for g in geoms)
        if d_first != d_last:
            flip = d_last < d_first
    neighbours = [n for n in G.successors(linkno) if n in streams_gdf.index]
    if flip is None and neighbours:
        geom = _first_geometry(streams_gdf, neighbours[0])
        d_first, d_last = first_point.distance(geom), last_point.distance(geom)
        if d_first != d_last:
            flip = d_first < d_last
    if flip is None:
        # Fallback to using elevation to determine upstream/downstream if no upstream or downstream nodes exist.
        flip = bool(dem[row1, col1] < dem[row2, col2])
    if flip:
        coords = np.flipud(coords)
        col1, row1, col2, row2 = col2, row2, col1, row1  # Swap the coordinates as well

    col2, row2 = col1, row1

    # Check if we have any upstream ids.
    if upstream_ids is None:
        upstream_ids = list(G.predecessors(linkno))
    last_elevation = np.inf
    starts = []
    for upstream_id in upstream_ids:
        us_row, us_col = nearest_stream_raster_pixel(streams, row1, col1, upstream_id)
        if 0 <= us_row < dem.shape[0] and 0 <= us_col < dem.shape[1] and dem[us_row, us_col] != nodata_value:
            starts.append(us_row * dem.shape[1] + us_col)
            last_elevation = min(last_elevation, dem[us_row, us_col])

    rows, cols = _line_cells(dem, streams, coords, linkno, inverse_transform, row1, col1, nodata_value)
    left = np.zeros(rows.size, np.bool_) if leveled is None else leveled[rows, cols]
    path = _burn_cells(dem, streams, rows, cols, left, local_min, last_elevation, nodata_value)
    if kept is not None:
        kept[rows[left], cols[left]] = True
    return (np.array(starts, np.int64), path) if path.size > 1 else None

def _first_geometry(streams_gdf: gpd.GeoDataFrame, stream_id):
    """A reach's geometry from a gdf indexed by stream id, taking the first row if the id is listed twice."""
    geom = streams_gdf.loc[stream_id, 'geometry']
    return geom.iloc[0] if isinstance(geom, pd.Series) else geom

@njit(cache=True, nogil=True)
def _burn_linestring(
    dem: np.ndarray,
    streams: np.ndarray,
    coords: np.ndarray,
    linkno: int,
    local_min: np.ndarray,
    gt: tuple,
    row1: int, 
    col1: int, 
    last_elevation: float,
    nodata_value: float) -> np.ndarray:
    """
    Lower the cells along ``coords`` (``gt`` is the inverse geotransform) to multiples of 0.5 m that never rise
    downstream, starting at (row1, col1) no higher than ``last_elevation``. Returns the cells, as flat indices, from
    upstream down.
    """
    rows, cols = _line_cells(dem, streams, coords, linkno, gt, row1, col1, nodata_value)
    return _burn_cells(dem, streams, rows, cols, np.zeros(rows.size, np.bool_), local_min, last_elevation,
                       nodata_value)


@njit(cache=True, nogil=True)
def _line_cells(
    dem: np.ndarray,
    streams: np.ndarray,
    coords: np.ndarray,
    linkno: int,
    gt: tuple,
    row1: int,
    col1: int,
    nodata_value: float) -> tuple[np.ndarray, np.ndarray]:
    """
    The cells along ``coords`` (``gt`` is the inverse geotransform), from (row1, col1) down: reach ``linkno``'s
    own where the stream raster has one next to a vertex, and the cells between vertices up to MAX_BURN_GAP_CELLS
    apart. Returns their rows and columns.
    """
    xs = coords[:, 0]
    ys = coords[:, 1]
    cols_ = np.floor(gt[0] + gt[1] * xs + gt[2] * ys).astype(np.int64)
    rows_ = np.floor(gt[3] + gt[4] * xs + gt[5] * ys).astype(np.int64)

    nrows, ncols = dem.shape

    rows = [row1]
    cols = [col1]
    last_row = row1
    last_col = col1

    for row, col in zip(rows_[1:], cols_[1:]):
        row, col = nearest_stream_raster_pixel(streams, row, col, linkno)

        if not (0 <= row < nrows and 0 <= col < ncols and dem[row, col] != nodata_value):
            continue

        if row != last_row or col != last_col:
            # Vertices more than a cell apart (a line drawn a cell's length in metres apart, where cells are
            # narrower than that one way) would skip the cells between them, and an unburned cell dams the channel.
            # Burn the cells between them too, the reach's own where the stream raster has one. A longer gap is
            # where the line left the DEM and came back, and is left alone.
            steps = max(abs(row - last_row), abs(col - last_col))
            if 1 < steps <= MAX_BURN_GAP_CELLS:
                for k in range(1, steps):
                    fr = last_row + (row - last_row) * k / steps
                    fc = last_col + (col - last_col) * k / steps
                    r = int(math.floor(fr + 0.5))
                    c = int(math.floor(fc + 0.5))
                    for rr in (int(math.floor(fr)), int(math.ceil(fr))):
                        for cc in (int(math.floor(fc)), int(math.ceil(fc))):
                            if 0 <= rr < nrows and 0 <= cc < ncols and streams[rr, cc] == linkno:
                                r, c = rr, cc
                    if 0 <= r < nrows and 0 <= c < ncols and dem[r, c] != nodata_value:
                        rows.append(r)
                        cols.append(c)
            rows.append(row)
            cols.append(col)
            last_row = row
            last_col = col
    return np.array(rows, np.int64), np.array(cols, np.int64)


@njit(cache=True, nogil=True)
def _burn_cells(
    dem: np.ndarray,
    streams: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
    keep: np.ndarray,
    local_min: np.ndarray,
    last_elevation: float,
    nodata_value: float) -> np.ndarray:
    """
    Lower a line's cells (``rows``, ``cols``, from upstream down) to multiples of 0.5 m that never rise downstream,
    starting no higher than ``last_elevation``, but for the ones to ``keep``: the line goes on down from those at
    their own level. Returns the cells, as flat indices, from upstream down.
    """
    nrows, ncols = dem.shape
    has_written = False

    # The trailing edge-of-DEM check below reads state from the last iteration of the loop.
    # A single-pixel reach never enters that loop, so these have to exist beforehand -- and
    # numba needs them bound on every path to compile the function at all. `has_written`
    # gates their use, so the seed values are never the ones acted on.
    row2 = rows[0]
    col2 = cols[0]

    # Every cell of the line is burned, even where it leaves the channel mask: an unburned gap dams the channel,
    # and filling turns the channel above it into a flat pool
    for k in range(rows.size - 1):
        row1, col1, row2, col2 = rows[k], cols[k], rows[k + 1], cols[k + 1]
        has_written = True

        upstream_elev = dem[row1, col1]
        downstream_elev = dem[row2, col2]

        if not keep[k]:
            if upstream_elev > last_elevation:
                upstream_elev = last_elevation

            min_non_stream_elev = local_min[row1, col1]
            if upstream_elev > min_non_stream_elev:
                upstream_elev = min_non_stream_elev

            if upstream_elev % 0.5 != 0:
                upstream_elev = math.floor(upstream_elev * 2) / 2

            dem[row1, col1] = upstream_elev

        if not keep[k + 1]:
            if downstream_elev > upstream_elev:
                downstream_elev = upstream_elev
            elif downstream_elev % 0.5 != 0:
                downstream_elev = math.floor(downstream_elev * 2) / 2

            dem[row2, col2] = downstream_elev
        last_elevation = downstream_elev
        if keep[k]:
            continue  # left as the DEM has it, and the ground around it too

        # Suppose that there are two+ neighbors, with an elevation <= new_elevation.
        # If one is in the stream raster but the other is not, let us bump the elevation of the other up
        minr = max(0, row1 - 1)
        maxr = min(nrows - 1, row1 + 1)
        minc = max(0, col1 - 1)
        maxc = min(ncols - 1, col1 + 1)
        rs, cs = np.nonzero(dem[minr:maxr+1, minc:maxc+1] <= upstream_elev)
        rs += minr
        cs += minc
        instream_neighbors = []
        outstream_neighbors = []
        for r, c in zip(rs, cs):
            if (r == row1 and c == col1) or dem[r, c] == nodata_value:
                continue
            if streams[r, c] > 0:
                instream_neighbors.append((r, c))
            else:
                outstream_neighbors.append((r, c))

        if instream_neighbors and outstream_neighbors:
            for r, c in outstream_neighbors:
                dem[r, c] = upstream_elev + 0.5

    # A reach that leaves the DEM, over its edge or into nodata (the ocean), ends 0.5 m lower, so the filled DEM routes
    # out there. Only there: each reach starts no higher than the reaches flowing into it, so lowering every reach's
    # end lowered the whole network below it 0.5 m at each junction, metres below a flat river's water downstream
    if has_written and dem[row2, col2] != nodata_value and (
            row2 == 0 or row2 == nrows - 1 or col2 == 0 or col2 == ncols - 1
            or (dem[row2 - 1:row2 + 2, col2 - 1:col2 + 2] == nodata_value).any()):
        dem[row2, col2] -= 0.5

    return rows * ncols + cols

# How far the D8 path may run outside a reach's cells before the reach is finished. A moved reach never leaves
# its D8 path, but a stream raster drawn from vector lines strays off it by a cell or so.
_MAX_OFF_REACH_STEPS = 2

@njit(cache=True)
def _trace_reach(cells: np.ndarray, flowdir: np.ndarray, nrows: int, ncols: int,
                 row_offsets: np.ndarray, col_offsets: np.ndarray,
                 max_off_reach_steps: int) -> tuple[int, int, int]:
    """
    Follow the D8 flow directions through one reach's cells (sorted flat indices).

    Returns ``(start_pixel, end_pixel, length)``: where the longest D8 path through the cells starts, its last
    cell in the reach, and how many cells it visits. Walking ``length`` cells down from ``start_pixel`` ends on
    ``end_pixel``.
    """
    n = cells.size

    # link[i] is the next cell of this reach that cells[i] drains into (-1 if the
    # path leaves the reach for good), and gap[i] is how many D8 steps that takes.
    # It is usually one step; more when the rasterized line strays off the D8 path.
    link = np.full(n, -1, np.int64)
    gap = np.zeros(n, np.int64)
    for i in range(n):
        pixel = cells[i]
        for step in range(1, max_off_reach_steps + 1):
            code = flowdir[pixel]
            if code < 0 or code > 255:
                break
            d_row = row_offsets[code]
            d_col = col_offsets[code]
            if d_row == 0 and d_col == 0:  # pit, outlet, or nodata
                break
            row = pixel // ncols + d_row
            col = pixel % ncols + d_col
            if row < 0 or row >= nrows or col < 0 or col >= ncols:
                break
            pixel = row * ncols + col
            j = np.searchsorted(cells, pixel)
            if j < n and cells[j] == pixel:
                link[i] = j
                gap[i] = step
                break

    # Every cell has at most one successor, so the reach forms a forest. Score each
    # cell by the path leading out of it and keep whichever root covers the most of
    # the reach: that is its head, and the far end of its path is its outlet.
    covered = np.zeros(n, np.int64)   # cells of this reach on the path out of i
    path_len = np.zeros(n, np.int64)  # D8 steps from i to the end of that path
    terminus = np.arange(n)
    state = np.zeros(n, np.uint8)     # 0 unvisited, 1 on the stack, 2 resolved
    stack = np.empty(n, np.int64)
    for root in range(n):
        if state[root] != 0:
            continue
        top = 0
        stack[0] = root
        state[root] = 1
        while top >= 0:
            i = stack[top]
            j = link[i]
            if j >= 0 and state[j] == 0:
                top += 1
                stack[top] = j
                state[j] = 1
                continue
            if j < 0 or state[j] == 1:  # end of the path, or a cycle across a flat
                link[i] = -1
                covered[i] = 1
                path_len[i] = 0
                terminus[i] = i
            else:
                covered[i] = 1 + covered[j]
                path_len[i] = gap[i] + path_len[j]
                terminus[i] = terminus[j]
            state[i] = 2
            top -= 1

    head = 0
    for i in range(1, n):
        if covered[i] > covered[head] or (covered[i] == covered[head] and path_len[i] > path_len[head]):
            head = i
    return cells[head], cells[terminus[head]], path_len[head] + 1

def _create_stream_info_table(
        stream_info_file: Path,
        streams_array: np.ndarray,
        flow_direction_file: str | Path) -> None:
    """
    Write the stream info table Curve2Flood's FLDPLN mapper reads.

    FLDPLN walks ``length`` cells down the D8 flow directions from each reach's ``start_pixel`` to find its
    stream cells, so each row must be a real D8 path that ends on ``end_pixel``. The table is therefore traced
    from the flow direction raster, not from the stream lines.
    """
    flow_direction_file = Path(flow_direction_file)
    if not flow_direction_file.exists():
        raise FileNotFoundError(
            f"Flow direction raster {flow_direction_file} does not exist. The FLDPLN stream info table is "
            "traced from it, so the hydrography has to be derived first (move_stream_network_to_thalweg)."
        )

    flowdir_ds: gdal.Dataset = gdal.Open(str(flow_direction_file))
    flowdir: np.ndarray = flowdir_ds.ReadAsArray()
    if flowdir.shape != streams_array.shape:
        raise ValueError(
            f"Flow direction raster {flow_direction_file} is {flowdir.shape} but the stream raster is "
            f"{streams_array.shape}. Both have to be on the DEM grid for the pixel indices to line up."
        )

    nrows, ncols = streams_array.shape
    linknos_flat = np.ascontiguousarray(streams_array).ravel()
    flowdir_flat = np.ascontiguousarray(flowdir).ravel().astype(np.int64, copy=False)
    row_offsets, col_offsets = _d8_offset_tables()

    # Group the stream pixels into one sorted block per link number.
    stream_pixels = np.flatnonzero(linknos_flat > 0)
    stream_pixels = stream_pixels[np.argsort(linknos_flat[stream_pixels], kind='stable')]
    linknos, block_starts = np.unique(linknos_flat[stream_pixels], return_index=True)
    block_ends = np.append(block_starts[1:], stream_pixels.size)

    output_table = []
    for linkno, block_start, block_end in zip(linknos, block_starts, block_ends):
        cells = np.sort(stream_pixels[block_start:block_end])
        start_pixel, end_pixel, length = _trace_reach(
            cells, flowdir_flat, nrows, ncols, row_offsets, col_offsets, _MAX_OFF_REACH_STEPS
        )
        output_table.append((int(start_pixel), int(end_pixel), int(length), int(linkno)))

    stream_info = pd.DataFrame(output_table,
                 columns=['start_pixel', 'end_pixel', 'length', 'stream_id'])
    stream_info_file.parent.mkdir(parents=True, exist_ok=True)
    if stream_info_file.suffix.lower() in {".parquet", ".pq"}:
        stream_info.to_parquet(stream_info_file, index=False)
    else:
        stream_info.to_csv(stream_info_file, index=False)


def make_channel_mask(channel_mask: np.ndarray, dem: str, water_mask: str, compression: str):
    """Write ``channel_mask`` to ``water_mask`` as a Byte raster on the DEM's grid."""
    dem_ds: gdal.Dataset = gdal.Open(dem)
    out_ds: gdal.Dataset = gdal.GetDriverByName('GTiff').Create(water_mask, dem_ds.RasterXSize, dem_ds.RasterYSize, 1, gdal.GDT_Byte, options=[f'COMPRESS={compression}'])
    out_ds.WriteArray(channel_mask)
    out_ds.SetGeoTransform(dem_ds.GetGeoTransform())
    out_ds.SetProjection(dem_ds.GetProjection())
    out_ds = None
