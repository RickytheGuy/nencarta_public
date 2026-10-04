import math
import shutil
import weakref
import tempfile
import warnings
from pathlib import Path

import tqdm
import numpy as np
import pandas as pd
import networkx as nx
from numba import njit
import geopandas as gpd
from osgeo import gdal, ogr, osr
from whitebox import WhiteboxTools
from scipy.ndimage import binary_dilation, distance_transform_edt, label, minimum_filter
from shapely.geometry import Point, LineString

from nencarta.logger import LOG
from nencarta.core.raster import Raster
from nencarta.core.vector import Vector
from nencarta.tasks.make_stream_geometry import _filter_streams_by_stream_order
from nencarta.tasks.stream_conflation import (
    ConflationSettings,
    conflate_to_flow_directions,
    D8_OFFSETS as _D8_OFFSETS,
    d8_offset_tables as _d8_offset_tables,
)
from nencarta.workspace import Workspace
from curve2flood import remove_cells_not_connected

def _in_memory_ogr_layer(gdf: gpd.GeoDataFrame, projection: str,
                         values: pd.Series | None = None) -> tuple[ogr.DataSource, ogr.Layer]:
    """
    Copy ``gdf``'s geometries, in order, into an in-memory OGR layer so GDAL can rasterize them, with ``values``
    in an integer field named "value" if given.

    The caller must keep the returned DataSource alive: OGR owns the layer, and letting the
    DataSource fall out of scope invalidates the layer while it is still in use.
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
    We load the lakes in the DEM's domain. We do not want to include lakes/reservoirs.

    The lakes are read through :class:`Vector` rather than ``ogr.Open`` so that every vector
    format nencarta accepts elsewhere works here too. ``ogr.Open`` cannot read a GeoParquet
    lakes file unless GDAL was built with the Parquet/Arrow driver; without it GDAL falls
    through to the ADBC driver and dies on a missing ``duckdb.dll``, which took out the whole
    burn/move step for the default (GeoParquet) lakes layer.
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
    if not workspace.configs.lakes:
        return None

    lakes_gdf = Vector(workspace.configs.lakes, not workspace.configs.parallel).to_geopandas(bbox_epsg_4326=dem_raster.epsg_4326_bbox)

    return lakes_gdf

def whitebox_callback(message: str) -> None:
    """
    Callback function for WhiteboxTools to log messa`ges.
    We only want to log errors and warnings, so we filter out other messages.
    """
    lowered = message.lower()
    if "warning" in lowered and not lowered.endswith("It appears that the input data is in"):
        LOG.warning(message)
    if "error" in lowered or "panic" in lowered:
        LOG.error(message) 

class _WhiteboxOutput:
    """Collects WhiteboxTools messages so a failure can report what the tool actually said."""

    def __init__(self, log: bool = True) -> None:
        self.lines: list[str] = []
        self.log = log

    def __call__(self, message: str) -> None:
        self.lines.append(message)
        if self.log:
            whitebox_callback(message)


def _private_whitebox(plugins: tuple[str, ...] = ()) -> WhiteboxTools:
    """
    Return a WhiteboxTools that runs its own copy of whitebox_tools (and of ``plugins``).

    whitebox_tools and its plugins read their settings from a settings.json next to the
    executable, and whitebox_tools rewrites that file whenever it is passed --compress_rasters or
    --max_procs, which the Python wrapper does on every tool run. When many processes share one
    install (a SLURM array, or a process pool), one reads the file while another is half way
    through rewriting it, and the tool panics with "Failed to parse config_file.json file".
    A private copy gives each run its own settings.json. Plugins are copied rather than
    symlinked because a symlinked executable resolves back to the shared install's settings.

    The copy is deleted when the returned object is garbage collected.
    """
    wbt = WhiteboxTools()
    shared_dir = Path(wbt.exe_path)
    private_dir = Path(tempfile.mkdtemp(prefix="nencarta_whitebox_"))
    weakref.finalize(wbt, shutil.rmtree, private_dir, True)

    shutil.copy2(shared_dir / wbt.exe_name, private_dir / wbt.exe_name)
    (private_dir / "plugins").mkdir()
    for plugin in plugins:
        exe = plugin + wbt.ext
        shutil.copy2(shared_dir / "plugins" / exe, private_dir / "plugins" / exe)
        shutil.copy2(shared_dir / "plugins" / f"{plugin}.json", private_dir / "plugins" / f"{plugin}.json")

    wbt.set_whitebox_dir(str(private_dir))
    return wbt


def _run_whitebox(tool, expected: Path, description: str, *args, **kwargs) -> int:
    """
    Run one WhiteboxTools tool and fail loudly if it crashed or did not produce ``expected``.

    The WhiteboxTools Python wrapper reads the tool's stdout until EOF and then returns 0
    unconditionally -- it never waits on the child or looks at its exit status. With verbose
    mode off the tool prints nothing at all, so a crash is indistinguishable from success and
    the only symptom is a missing output file. whitebox_tools.exe does crash: when processes
    share one install it panics (exit code 101) on a minority of runs (see _private_whitebox),
    and without the child's real exit status there is nothing to tell that apart from a bad input. Tools that rewrite their input
    in place leave ``expected`` behind even when they crash, so the exit status is always checked.

    The wrapper also only hands the tool's output to the callback in verbose mode, so with
    verbose off a panic message is thrown away. Verbose is forced on for the call so the output
    can be reported on failure; it is only logged if the caller had verbose on.
    """
    import whitebox.whitebox_tools as whitebox_module

    wbt = tool.__self__
    was_verbose = wbt.verbose
    output = _WhiteboxOutput(log=was_verbose)
    launched = []
    original_popen = whitebox_module.Popen

    def recording_popen(*popen_args, **popen_kwargs):
        process = original_popen(*popen_args, **popen_kwargs)
        launched.append(process)
        return process

    whitebox_module.Popen = recording_popen
    # Set the attribute directly: set_verbose_mode() also rewrites the shared settings.json
    wbt.verbose = True
    try:
        code = tool(*args, callback=output, **kwargs)
    finally:
        whitebox_module.Popen = original_popen
        wbt.verbose = was_verbose

    exit_codes = []
    for process in launched:
        try:
            exit_codes.append(process.wait(timeout=60))
        except Exception as exc:
            exit_codes.append(f"<{type(exc).__name__}>")

    if expected.exists() and all(exit_code == 0 for exit_code in exit_codes):
        return code

    reported = "\n".join(output.lines[-15:]) or "<the tool printed nothing>"
    panicked = " (101 is a Rust panic)" if 101 in exit_codes else ""
    error = RuntimeError if expected.exists() else FileNotFoundError
    raise error(
        f"{description} was not created successfully: {expected}. "
        f"whitebox_tools exited with {exit_codes}{panicked}. Output:\n{reported}"
    )

# The only GeoTIFF compressions the WhiteboxTools decoder can read
_WHITEBOX_TIFF_COMPRESSIONS = {"NONE", "PACKBITS", "LZW", "DEFLATE"}

def _whitebox_readable_dem(dem_path: Path, scratch_path: Path) -> Path:
    """
    Return a path to ``dem_path`` that WhiteboxTools can read.

    WhiteboxTools has its own GeoTIFF decoder that only supports PACKBITS, LZW, and DEFLATE,
    and panics (exit code 101) on anything else. The burned DEM is written with
    ``configs.compression``, so a ``compression: ZSTD`` run (or a VRT input) would crash the
    fill step. In that case the DEM is copied to a DEFLATE GeoTIFF at ``scratch_path``.
    """
    ds: gdal.Dataset = gdal.Open(str(dem_path))
    compression = (ds.GetMetadataItem("COMPRESSION", "IMAGE_STRUCTURE") or "NONE").upper()
    is_geotiff = ds.GetDriver().ShortName == "GTiff"
    ds = None
    if is_geotiff and compression in _WHITEBOX_TIFF_COMPRESSIONS:
        return dem_path

    gdal.Translate(str(scratch_path), str(dem_path), format="GTiff", creationOptions=["COMPRESS=DEFLATE"])
    return scratch_path


def derive_flow_directions_using_whitebox(workspace: Workspace, dem_path: Path) -> None:
    """
    Fill the depressions in ``dem_path`` and derive its D8 flow directions with WhiteboxTools, writing
    workspace.filled_dem and workspace.flowdir. The stream network is moved onto these flow directions, and
    FLDPLN floods along them.
    """
    wbt = _private_whitebox()
    wbt.set_compress_rasters(True)
    wbt.set_verbose_mode(LOG.level <= 20)  # INFO or lower
    wbt.set_max_procs(1)
    workspace.dem_updated_folder.mkdir(parents=True, exist_ok=True)
    workspace.Flow_Direction_Folder.mkdir(parents=True, exist_ok=True)

    # Though fill_depressions is more efficient thanfill_depressions_wang_and_liu,
    # fill_depressions_wang_and_liu is more stable and crashes less
    dem_path = Path(dem_path)
    whitebox_dem = _whitebox_readable_dem(
        dem_path,
        workspace.dem_updated_folder / f"{dem_path.stem}_whitebox.tif"
    )
    try:
        _run_whitebox(wbt.fill_depressions_wang_and_liu, workspace.filled_dem,
                    f"Filled DEM from {dem_path}",
                    str(whitebox_dem), str(workspace.filled_dem))
    finally:
        if whitebox_dem != dem_path:
            whitebox_dem.unlink(missing_ok=True)
    _run_whitebox(wbt.d8_pointer, workspace.flowdir, "Flow direction file",
                  str(workspace.filled_dem), str(workspace.flowdir))

def burn_streams_and_move_streams(workspace: Workspace) -> Path:
    """
    This function burns streams into the DEM and/or moves the streams onto the DEM's D8 flow paths.
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

    if should_burn_streams:
        workspace.dem_updated_folder.mkdir(parents=True, exist_ok=True)

        lakes = load_lake_array(workspace, assigned_dem)

        source_gdf = Vector(workspace.DEM_StrmShp, not configs.parallel).to_geopandas().to_crs(assigned_dem.projection)
        channel_mask, dem_for_conflation = smooth_and_burn_dem(
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

        channel_mask = Raster(workspace.bathy_water_mask).read_array()
        source_gdf = Vector(workspace.DEM_StrmShp, not configs.parallel).to_geopandas().to_crs(assigned_dem.projection)

    if should_move_streams:
        derive_flow_directions_using_whitebox(workspace, dem_for_conflation_path)
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
        _write_stream_raster(workspace.new_stream_raster, final_streams, assigned_dem)

        channel_mask |= (final_streams > 0)
        channel_mask = remove_cells_not_connected(channel_mask, final_streams)

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
    Move the source streams onto workspace.flowdir's D8 paths (see stream_conflation) and return them with the
    stream raster they paint.

    Only DEM channels draining at least new_strm_threshold_km2 are matched onto. Lakes are cut out for FLDPLN,
    which isn't meant to map them: no stream runs inside one, and a reach flowing into one ends at its shore. The
    other mappers keep the old rule, dropping the reaches that touch a lake unless they carry a river through it.
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
    area_col = None
    if configs.area_km2_field and configs.area_km2_field in source_gdf.columns:
        area_col = configs.area_km2_field
    elif configs.area_m2_field and configs.area_m2_field in source_gdf.columns:
        area_col = "__area_km2"
        source_gdf = source_gdf.assign(__area_km2=pd.to_numeric(source_gdf[configs.area_m2_field], errors="coerce") / 1e6)
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
        # A reach whose source line can't be merged into one line is drawn in several pieces. Whether it is has to
        # be decided after a line_merge, not read off the stored geometry type: a GeoPackage layer stores every
        # LineString as a MultiLineString, so on the N14W89 domain the stored type flagged all 665 reaches when
        # only 165 are really split.
        merged = source_gdf.geometry.line_merge()
        multipart = set(source_gdf.loc[(merged.geom_type == 'MultiLineString').to_numpy(), id_col])
        streams, raster = _drop_reaches(streams, raster, multipart, id_col, ds_col)

    if configs.StrmOrder_Field and (configs.StrmOrder_Lower is not None or configs.StrmOrder_Upper is not None) and not configs.mapper.is_curve2flood_fldpln_mapper():
        kept = _filter_streams_by_stream_order(streams.copy(), configs.StrmOrder_Field, configs.StrmOrder_Lower, configs.StrmOrder_Upper)
        streams, raster = _drop_reaches(streams, raster, set(streams[id_col]) - set(kept[id_col]), id_col, ds_col)

    return streams, raster

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

    if workspace.configs.use_dem_derived_channel_mask:
        channel_mask = (dem % 0.5 == 0)
    else:
        channel_mask = Raster(workspace.bathy_water_mask).read_array()

    channel_mask = burn_streams_into_dem(channel_mask, dem, streams, dem_ds, source_gdf, lakes, id_col, ds_col, nodata_value)
    dem = smooth_burned_dem(dem, channel_mask, streams, pbar=False)

    dem[dem < -1000] = nodata_value # Remove any DEM values that are less than -1000 m, since these are likely to be erroneous and will cause problems with the floodplain mapping.
    output_ds = gdal.GetDriverByName('GTiff').Create(workspace.fixed_dem, dem_ds.RasterXSize, dem_ds.RasterYSize, 1, gdal.GDT_Float32, options=[f'COMPRESS={workspace.configs.compression}'])
    output_ds.WriteArray(dem)
    output_ds.SetGeoTransform(dem_ds.GetGeoTransform())
    output_ds.SetProjection(dem_ds.GetProjection())
    output_ds.GetRasterBand(1).SetNoDataValue(nodata_value)
    output_ds = None

    return channel_mask, dem

def build_mask_graph(dem: np.ndarray, mask: np.ndarray) -> nx.Graph:
    G = nx.Graph()

    rows, cols = np.nonzero(mask)

    G.add_nodes_from(
        ((r, c), {"elevation": dem[r, c]})
        for r, c in zip(rows, cols)
    )

    # Horizontal
    rr, cc = np.nonzero(mask[:, :-1] & mask[:, 1:])
    G.add_edges_from(zip(zip(rr, cc), zip(rr, cc + 1)))

    # Vertical
    rr, cc = np.nonzero(mask[:-1, :] & mask[1:, :])
    G.add_edges_from(zip(zip(rr, cc), zip(rr + 1, cc)))

    # Diagonal
    rr, cc = np.nonzero(mask[:-1, :-1] & mask[1:, 1:])
    G.add_edges_from(zip(zip(rr, cc), zip(rr + 1, cc + 1)))

    return G


def elevation_components_nx(G: nx.Graph, elevations: dict) -> list[set]:
    visited = set()
    components = []

    # Traverse graph, finding connected components of nodes with the same elevation. Each component is a set of (row, col) tuples.
    for node in G.nodes:
        if node in visited:
            continue

        elev = elevations[node]
        stack = [node]
        component = set()

        while stack:
            current = stack.pop()
            if current in visited:
                continue

            visited.add(current)
            component.add(current)

            for neighbor in G.neighbors(current):
                if neighbor not in visited and elevations[neighbor] == elev:
                        stack.append(neighbor)

        components.append(component)

    return components

def smooth_burned_dem(dem: np.ndarray, mask: np.ndarray = None, streams: np.ndarray = None, pbar: bool = True, max_difference: float = 0.5) -> np.ndarray:
    dem = dem.astype(np.float32, copy=False)
    if mask is None:
        mask = (dem % max_difference) == 0

    banks = binary_dilation(mask, structure=np.ones((3, 3), dtype=int)).astype(np.uint8, copy=False)
    banks[mask] = False  # Banks are the cells adjacent to the mask, but not in the mask

    # Get a distance raster for the mask, where the distance is highest in the center, but smallest at the edges
    distance_to_banks = distance_transform_edt(mask)

    G = build_mask_graph(dem, mask)
    elevations = nx.get_node_attributes(G, "elevation")
    components = elevation_components_nx(G, elevations)
    final_node_elevations = {}
    for component in tqdm.tqdm(components, desc="Processing elevation components", disable=not pbar):
        if len(component) <= 1:
            continue

        # 2: Find all nodes that are on the outside edge of the component (a node that outside component but has a neighbor inside the component, but within the mask)
        boundary_nodes = set()
        max_distance_to_banks = 0
        min_distance_to_banks = float('inf')
        for node in component:
            for neighbor in G.neighbors(node):
                if neighbor not in component:
                    boundary_nodes.add(neighbor)
                    if distance_to_banks[neighbor] > max_distance_to_banks:
                        max_distance_to_banks = distance_to_banks[neighbor]
                    if distance_to_banks[neighbor] < min_distance_to_banks:
                        min_distance_to_banks = distance_to_banks[neighbor]

        # 3: Starting at the boundary, encroach into the component, and for each node, compute a new elevation
        # as the number of component nodes encroached over the total nodes times the difference between the boundary elevation and the starting elevation.
        if not boundary_nodes:
            continue

        modified_node_elevations = {}
        source_elevation = elevations[list(component)[0]]

        # 3.1 Find all nodes in the component that are attatched to a boundary node
        visited = set()
        global_index = 0
        while boundary_nodes:
            new_boundary = set()
            for node in boundary_nodes:
                for neighbor in G.neighbors(node):
                    if neighbor in component:
                        if neighbor in visited:
                            continue
                        if global_index == 0 and abs(elevations[node] - source_elevation) > max_difference:
                            continue
                        new_boundary.add(neighbor)
                        visited.add(neighbor)
                        if global_index == 0:
                            modified_node_elevations[neighbor] = (elevations[node], global_index)
                        else:
                            modified_node_elevations[neighbor] = (modified_node_elevations[node][0], global_index)

            global_index += 1
            boundary_nodes = new_boundary

        # Calculate the max index for each boundary elevation
        boundary_elevation_to_max_index = {}
        for node, (boundary_elevation, distance) in modified_node_elevations.items():
            if boundary_elevation not in boundary_elevation_to_max_index:
                boundary_elevation_to_max_index[boundary_elevation] = distance
            else:
                boundary_elevation_to_max_index[boundary_elevation] = max(boundary_elevation_to_max_index[boundary_elevation], distance)

        for node, (boundary_elevation, distance) in modified_node_elevations.items():
            local_index = boundary_elevation_to_max_index[boundary_elevation]
            if local_index == 0:
                continue
            else:
                # Interpolation, which does not exactly pass through neither the boundary nor the source elevation
                # It will bring us close to those values, which allows for a smoother transition rather than
                # two of the same elevations back to back across components or within a component across internal meeting boundary.
                new_elevation = boundary_elevation + ((source_elevation - boundary_elevation) / 2) + (source_elevation - boundary_elevation) * ((distance+1) / (local_index+2) / 2)
            
                # Lower elevation up to 0.1 m to help stream be in the center of the channel. ARC allows flat water detection up to 0.1 m
                if max_distance_to_banks != min_distance_to_banks:
                    new_elevation -= (0.1 * (max(min(distance_to_banks[node], max_distance_to_banks), min_distance_to_banks) - min_distance_to_banks) / (max_distance_to_banks - min_distance_to_banks))

            # Ensure that the new elevation does not exceed the banks, by raising the banks a bit
            dem_mask = (dem[node[0]-1:node[0]+2, node[1]-1:node[1]+2] <= new_elevation)
            rows, cols = np.nonzero(
                dem_mask & banks[node[0]-1:node[0]+2, node[1]-1:node[1]+2]
            )
            if len(rows) > 0:
                dem[rows + node[0]-1, cols + node[1]-1] = new_elevation + 0.5

            final_node_elevations[node] = new_elevation

    for node, new_elevation in final_node_elevations.items():
        dem[node] = new_elevation

    return dem

def _get_stream_raster(dem_ds: gdal.Dataset, streams_gdf: gpd.GeoDataFrame, id_col: str = 'LINKNO',
                       order_col: str | None = None, area_col: str | None = None) -> np.ndarray:
    """
    Rasterize the stream ids onto the DEM's grid.

    Where reaches meet, a cell can carry only one of their ids, the last one drawn, and the burn lowers it only
    as far as that reach goes. Drawn in table order, the cell where a tributary joins a river often went to the
    tributary and was left at the tributary's level: a dam across the river, behind which fill turns the burned
    channel into a flat pool. On the N38W085 tile one such cell flattened 30 km of the Ohio, and the D8 flow
    across that flat ran in straight lines instead of down the river. The larger reaches, by stream order and then
    drainage area, are therefore drawn last, so that a river keeps the cells it shares with its tributaries.
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
        min_feature_size: int = 5
):
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

    # Burn upstream reaches first. A reach starts no higher than the reaches flowing into it (see burn_linestring),
    # which only holds if they were burned before it; in table order a river could be burned before its
    # tributaries had been, leaving a step up where they meet.
    try:
        burn_order = {stream_id: i for i, stream_id in enumerate(nx.topological_sort(G))}
        streams_gdf = streams_gdf.iloc[np.argsort(streams_gdf[id_col].map(burn_order).to_numpy(), kind="stable")]
    except nx.NetworkXUnfeasible:
        LOG.warning("The source stream network has a cycle, so its reaches are burned in table order.")
    streams_gdf = streams_gdf.set_index(id_col)

    masked_dem = np.where(channel_border, dem, np.inf)
    local_min = minimum_filter(masked_dem, size=3, mode="nearest")

    # Traverse each segment. Identify which end is up and downstream.
    # Then, take the upstream value. If it is not a multiple of 0.5, lower it to so.
    # Go downstream. Each cell is a multiple of 0.5, no higher than upstream elevation or the current cell elevation.
    for row in streams_gdf.itertuples():
        geom = row.geometry
        stream_id = row.Index
        if geom.geom_type == "MultiLineString":
            for line in geom.geoms:
                burn_linestring(dem, streams, dem_ds, line, stream_id, G, streams_gdf, local_min, nodata_value)
        elif geom.geom_type == "LineString":
            burn_linestring(dem, streams, dem_ds, geom, stream_id, G, streams_gdf, local_min, nodata_value)
        else:
            raise ValueError(f"Unsupported geometry type: {geom.geom_type}")
        
    return channel_mask

@njit(cache=True, nogil=True)
def nearest_stream_raster_pixel(
    streams_array: np.ndarray,
    start_row: int,
    start_col: int,
    linkno: int
):
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
        nodata_value: float):
    if linkno not in G:
        return  # Skip if the linkno is not in the graph
    
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
            return  # No valid point found within DEM bounds
        
    if not (0 <= row2 < dem.shape[0] and 0 <= col2 < dem.shape[1] and dem[math.floor(row2), math.floor(col2)] != nodata_value):
        # Traverse the linestring in reverse to find the last point that is within the DEM bounds
        for x, y in np.flipud(coords):
            col, row = gdal.ApplyGeoTransform(inverse_transform, x, y)
            if 0 <= row < dem.shape[0] and 0 <= col < dem.shape[1] and dem[math.floor(row), math.floor(col)] != nodata_value:
                col2, row2 = col, row
                break
        else:
            return  # No valid point found within DEM bounds
        
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
    upstream_ids = list(G.predecessors(linkno))
    last_elevation = np.inf
    for upstream_id in upstream_ids:
        us_row, us_col = nearest_stream_raster_pixel(streams, row1, col1, upstream_id)
        if 0 <= us_row < dem.shape[0] and 0 <= us_col < dem.shape[1] and dem[us_row, us_col] != nodata_value and dem[us_row, us_col] < last_elevation:
            last_elevation = dem[us_row, us_col]

    _burn_linestring(dem, streams, coords, linkno, local_min, inverse_transform, row1, col1, last_elevation, nodata_value)

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
    nodata_value: float) -> None:
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
            rows.append(row)
            cols.append(col)
            last_row = row
            last_col = col

    has_written = False

    # The trailing edge-of-DEM check below reads state from the last iteration of the loop.
    # A single-pixel reach never enters that loop, so these have to exist beforehand -- and
    # numba needs them bound on every path to compile the function at all. `has_written`
    # gates both uses, so the seed values are never the ones acted on.
    row2 = row1
    col2 = col1
    outstream_neighbor_count = 0

    # Every cell of the line is burned. This used to skip the stretches where the line left the channel mask and
    # came back into the same piece of it, on the idea that the water body carries the flow between them. It
    # can't: the burned channel sits up to 0.5 m below the water's surface, so water leaving it has to rise to the
    # surface to get past the gap, and fill turns the channel above the gap into a flat pool. Nor were the gaps
    # mostly real detours: away from land-cover water the mask is the cleaned stream raster, which differs from the
    # stream raster burned here by a cell in places, and each such one-cell "detour" dammed the river. On the
    # N38W085 tile the skips left half of the order 6+ river cells in pools that fill flattened.
    for row1, col1, row2, col2 in zip(rows[:-1], cols[:-1], rows[1:], cols[1:]):
        has_written = True

        upstream_elev = dem[row1, col1]
        downstream_elev = dem[row2, col2]

        if upstream_elev > last_elevation:
            upstream_elev = last_elevation

        min_non_stream_elev = local_min[row1, col1]
        if upstream_elev > min_non_stream_elev:
            upstream_elev = min_non_stream_elev

        if upstream_elev % 0.5 != 0:
            upstream_elev = math.floor(upstream_elev * 2) / 2

        dem[row1, col1] = upstream_elev

        if downstream_elev > upstream_elev:
            downstream_elev = upstream_elev
        elif downstream_elev % 0.5 != 0:
            downstream_elev = math.floor(downstream_elev * 2) / 2

        dem[row2, col2] = downstream_elev
        last_elevation = downstream_elev

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

        outstream_neighbor_count = len(outstream_neighbors)

        if instream_neighbors and outstream_neighbors:
            for r, c in outstream_neighbors:
                dem[r, c] = upstream_elev + 0.5

    # One more thing: check if the last (row2, col2) is on the border of the dem. If so, drop by 0.5 (helps filled dem step route out of the DEM)
    if has_written and (row2 == 0 or row2 == nrows - 1 or col2 == 0 or col2 == ncols - 1) and dem[row2, col2] != nodata_value:
        dem[row2, col2] -= 0.5
    # Same check, but for if we are on the edge of nodata
    elif has_written and outstream_neighbor_count == 0 and dem[row2, col2] != nodata_value:
        dem[row2, col2] -= 0.5

# How far the D8 path may run outside a reach's rasterized cells before we call the
# reach finished. The conflated stream raster is painted along D8 paths, so its reaches
# never leave them; a stream raster rasterized from vector lines does, by a cell or so
# (over 95% of the gaps measured on real tiles are a single cell), and a short tolerance
# keeps such a reach whole while still stopping a badly placed one from running away down
# the network.
_MAX_OFF_REACH_STEPS = 2

@njit(cache=True)
def _trace_reach(cells: np.ndarray, flowdir: np.ndarray, nrows: int, ncols: int,
                 row_offsets: np.ndarray, col_offsets: np.ndarray,
                 max_off_reach_steps: int) -> tuple[int, int, int]:
    """
    Follow the D8 flow direction through one reach's cells.

    ``cells`` is the sorted flat index of every pixel carrying this reach's link
    number. Returns ``(start_pixel, end_pixel, length)``: the pixel the longest D8
    path through those cells begins at, the last cell of the reach on that path, and
    the number of pixels the path visits, inclusive. Walking ``length`` pixels
    downstream from ``start_pixel`` therefore lands exactly on ``end_pixel``.
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
    Write the reach table that Curve2Flood's FLDPLN spreader reads.

    For each reach the spreader takes ``start_pixel`` and walks ``length`` pixels down
    the D8 flow direction raster, treating what it visits as that reach's stream
    pixels, then traces from ``end_pixel`` to exclude everything further downstream.
    Every row therefore has to describe a real D8 path. Taking ``length`` to be the
    number of rasterized cells carrying the link number, and the two endpoints from
    the ends of the stream vector, does not: rasterizing a line and tracing D8 across
    a filled DEM disagree about which cells belong to the reach, so the walk drifts
    off the reach and stops somewhere other than ``end_pixel``. The table is read off
    the flow direction raster instead, which is the same thing the spreader walks.
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
    dem_ds: gdal.Dataset = gdal.Open(dem)
    out_ds: gdal.Dataset = gdal.GetDriverByName('GTiff').Create(water_mask, dem_ds.RasterXSize, dem_ds.RasterYSize, 1, gdal.GDT_Byte, options=[f'COMPRESS={compression}'])
    out_ds.WriteArray(channel_mask)
    out_ds.SetGeoTransform(dem_ds.GetGeoTransform())
    out_ds.SetProjection(dem_ds.GetProjection())
    out_ds = None
