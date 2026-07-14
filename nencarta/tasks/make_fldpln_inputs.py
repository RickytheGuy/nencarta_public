import os
import math
import warnings
from enum import Enum
from pathlib import Path
from functools import partial
from dataclasses import dataclass
from itertools import combinations
from collections import defaultdict

import tqdm
import numpy as np
import pandas as pd
import networkx as nx
import geopandas as gpd
from numba import njit, prange
from osgeo import gdal, ogr
from whitebox import WhiteboxTools
from scipy.interpolate import griddata
from shapely import reverse, line_merge, prepare
from scipy.ndimage import binary_dilation, distance_transform_edt, label, minimum_filter
from shapely.geometry import box, Point, LineString, Polygon, MultiLineString, GeometryCollection

from nencarta.logger import LOG
from nencarta.core.raster import Raster
from nencarta.core.vector import Vector
from nencarta.workspace import Workspace
from nencarta.tasks.configs import define_arc_configs_for_fldpln_bathymetry
from arc import Arc
from curve2flood import remove_cells_not_connected

wbt = WhiteboxTools()
wbt.set_compress_rasters(True)

def make_fldpln_inputs(workspace: Workspace) -> Path:
    if workspace.stream_info_file.exists() and \
        workspace.filled_dem.exists() and \
        workspace.new_stream_raster.exists() and \
        workspace.new_StrmShp_matched.exists() and \
        workspace.flowdir.exists() and \
        workspace.fldpln_bathymetry.exists() and \
        not workspace.configs.process_stream_network:
        LOG.info(f"{workspace.stream_info_file} already exists and we aren't making it again...")
        workspace.assigned_dem = workspace.fldpln_bathymetry
        workspace.DEM_StrmShp = workspace.new_StrmShp_matched
        workspace.STRM_File_Clean = workspace.new_stream_raster
        workspace.configs.disable_bathymetry = True
        return workspace.stream_info_file
    
    if not workspace.DEM_StrmShp.exists() and not workspace.configs.raise_errors_if_nothing_in_domain:
        return None
    
    wbt.set_verbose_mode(LOG.level <= 20)  # INFO or lower
    
    workspace.dem_updated_folder.mkdir(parents=True, exist_ok=True)
    workspace.Flow_Direction_Folder.mkdir(parents=True, exist_ok=True)

    dem_raster = Raster(workspace.assigned_dem)
    if workspace.configs.lakes:
        lakes_ds: gdal.Dataset = gdal.GetDriverByName('MEM').Create('', dem_raster.shape[1], dem_raster.shape[0], 1, gdal.GDT_Byte)
        lakes_ds.SetGeoTransform(dem_raster.geotransform)
        lakes_ds.SetProjection(dem_raster.projection)
        ds: gdal.Dataset = ogr.Open(workspace.configs.lakes)
        lakes_layer: ogr.Layer = ds.GetLayer()
        bbox = dem_raster.bbox
        lakes_layer.SetSpatialFilterRect(*bbox)
        gdal.RasterizeLayer(lakes_ds, [1], lakes_layer, burn_values=[1])
        lakes_ds.FlushCache()
        lakes = lakes_ds.ReadAsArray().astype(np.bool_, copy=False)
    else:
        lakes = None

    source_gdf = Vector(workspace.DEM_StrmShp, not workspace.configs.parallel).to_geopandas()
    channel_mask, fixed_dem = smooth_and_burn_dem(
        workspace, 
        source_gdf, 
        lakes,
        workspace.configs.streamflow_source.upstream_id, 
        workspace.configs.streamflow_source.downstream_id
    )

    wbt.fill_depressions(str(workspace.fixed_dem), str(workspace.filled_dem))
    if not workspace.filled_dem.exists():
        raise FileNotFoundError(f"Filled DEM file {workspace.filled_dem} was not created successfully.")
    wbt.d8_pointer(str(workspace.filled_dem), str(workspace.flowdir))
    if not workspace.flowdir.exists():
        raise FileNotFoundError(f"Flow direction file {workspace.flowdir} was not created successfully.")
    wbt.d8_flow_accumulation(str(workspace.flowdir), str(workspace.flowacc), pntr=True, out_type='catchment area')
    if not workspace.flowacc.exists():
        raise FileNotFoundError(f"Flow accumulation file {workspace.flowacc} was not created successfully.")

    # # If DEM units are not in km2, convert the threshold to the DEM units. This is important for the stream extraction step, which uses the flow accumulation raster to determine where streams are.
    threshold = workspace.configs.new_strm_threshold_km2
    threshold_native = threshold * dem_raster.native_cell_area / dem_raster.cell_area_km2
    buffer_distance = (50/1000) * np.sqrt(dem_raster.native_cell_area / dem_raster.cell_area_km2) # 50 m buffer distance in DEM units
    snap_distance = (0.1/1000) * np.sqrt(dem_raster.native_cell_area / dem_raster.cell_area_km2) # 0.1 m snap distance in DEM units

    # Remove the vector rasters, since whitebox will not make some of them if they exist
    for file in workspace.new_StrmShp.parent.glob("*"):
        if file.stem.startswith(workspace.new_StrmShp.stem):
            file.unlink()

    wbt.extract_streams(str(workspace.flowacc), str(workspace.whitebox_stream_raster), threshold=threshold_native, zero_background=True)
    if not workspace.whitebox_stream_raster.exists():
        raise FileNotFoundError(f"Whitebox stream raster file {workspace.whitebox_stream_raster} was not created successfully. The threshold used was {threshold_native} in DEM units, which is equivalent to {threshold} km2.")
            
    # Remove in the newly created stream raster streams where the dem == 0
    streams_ds: gdal.Dataset = gdal.Open(str(workspace.whitebox_stream_raster), gdal.GA_Update)
    streams = streams_ds.ReadAsArray()
    streams[fixed_dem == 0] = 0
    streams_ds.WriteArray(streams)
    streams_ds = None

    wbt.stream_link_identifier(str(workspace.flowdir), str(workspace.whitebox_stream_raster), str(workspace.whitebox_stream_raster), zero_background=True)
    wbt.raster_streams_to_vector(str(workspace.whitebox_stream_raster), str(workspace.flowdir), str(workspace.new_StrmShp))
    if not workspace.new_StrmShp.exists():
        raise FileNotFoundError(f"New stream shapefile {workspace.new_StrmShp} was not created successfully.")

    # Whitebox does not insert the projection into the shapefile, so we need to do that here.
    prj_file = workspace.new_StrmShp.with_suffix('.prj')
    with open(prj_file, 'w') as f:
        f.write(dem_raster.projection)

    if lakes is not None:
        lakes_gdf = Vector(workspace.configs.lakes, not workspace.configs.parallel).to_geopandas(bbox_epsg_4326=dem_raster.epsg_4326_bbox)
        wtbx_gdf = Vector(workspace.new_StrmShp).to_geopandas()
        wtbx_gdf = wtbx_gdf[~wtbx_gdf.geometry.intersects(lakes_gdf.union_all())]  # Remove streams that intersect lakes
        wtbx_gdf.to_file(workspace.new_StrmShp)

    wbt.repair_stream_vector_topology(str(workspace.new_StrmShp), str(workspace.new_StrmShp), dist=snap_distance)
    wbt.vector_stream_network_analysis(str(workspace.new_StrmShp), str(workspace.filled_dem), str(workspace.new_StrmShp), snap=snap_distance)

    streams_gdf = _conflate_streams(
        source_gdf=source_gdf, 
        streams_vector=workspace.new_StrmShp, 
        buffer_distance=buffer_distance,
        dem_proj=dem_raster.projection, 
        dem_bbox=dem_raster.bbox,
        source_id_col=workspace.configs.streamflow_source.upstream_id,
        source_ds_col=workspace.configs.streamflow_source.downstream_id,
    )
    if streams_gdf.empty:
        if workspace.configs.raise_errors_if_nothing_in_domain:
            raise ValueError("No stream geometries remain after conflation.")
        else:
            workspace.DEM_StrmShp = workspace.new_StrmShp_matched
            return None
        
    Vector.save_any_geom(streams_gdf, workspace.new_StrmShp_matched)

    _rasterize_streams(str(workspace.new_stream_raster), str(workspace.fixed_dem), str(workspace.new_StrmShp_matched), attribute=workspace.configs.streamflow_source.upstream_id)

    final_streams = gdal.Open(str(workspace.new_stream_raster)).ReadAsArray()
    channel_mask |= (final_streams > 0)
    channel_mask = remove_cells_not_connected(channel_mask, final_streams)

    if not workspace.configs.disable_bathymetry:
        config = define_arc_configs_for_fldpln_bathymetry(workspace)
        Arc(str(config), quiet=workspace.configs.quiet).run()
        if workspace.ARC_BathyFile.exists():
            for file in [workspace.flowdir, workspace.flowacc, workspace.whitebox_stream_raster, workspace.new_StrmShp]:
                file.unlink()

            bathy_raster = Raster(workspace.ARC_BathyFile)
            bathy = bathy_raster.read_array()

            filled = interpolate_bathymetry(bathy, channel_mask)

            smoothed = smooth_downhill(
                filled,
                channel_mask,
                alpha=4.0,
                iterations=8,
            )
            smoothed[~channel_mask] = Raster(workspace.fixed_dem).read_array()[~channel_mask]

            Raster.write_array_using_reference(
                smoothed,
                bathy_raster,
                workspace.fldpln_bathymetry
            )

            wbt.fill_depressions(str(workspace.fldpln_bathymetry), str(workspace.filled_dem))
            if not workspace.filled_dem.exists():
                raise FileNotFoundError(f"Filled DEM file {workspace.filled_dem} was not created successfully.")
            wbt.d8_pointer(str(workspace.filled_dem), str(workspace.flowdir))
            if not workspace.flowdir.exists():
                raise FileNotFoundError(f"Flow direction file {workspace.flowdir} was not created successfully.")
            wbt.d8_flow_accumulation(str(workspace.flowdir), str(workspace.flowacc), pntr=True, out_type='catchment area')
            if not workspace.flowacc.exists():
                raise FileNotFoundError(f"Flow accumulation file {workspace.flowacc} was not created successfully.")

            # Remove the vector rasters, since whitebox will not make some of them if they exist
            for file in workspace.new_StrmShp.parent.glob("*"):
                if file.stem.startswith(workspace.new_StrmShp.stem):
                    file.unlink()

            wbt.extract_streams(str(workspace.flowacc), str(workspace.whitebox_stream_raster), threshold=threshold_native, zero_background=True)
            if not workspace.whitebox_stream_raster.exists():
                raise FileNotFoundError(f"Whitebox stream raster file {workspace.whitebox_stream_raster} was not created successfully. The threshold used was {threshold_native} in DEM units, which is equivalent to {threshold} km2.")
            wbt.stream_link_identifier(str(workspace.flowdir), str(workspace.whitebox_stream_raster), str(workspace.whitebox_stream_raster), zero_background=True)
            wbt.raster_streams_to_vector(str(workspace.whitebox_stream_raster), str(workspace.flowdir), str(workspace.new_StrmShp))
            if not workspace.new_StrmShp.exists():
                if workspace.configs.raise_errors_if_nothing_in_domain:
                    raise FileNotFoundError(f"New stream shapefile {workspace.new_StrmShp} was not created successfully.")
                else:
                    if workspace.new_StrmShp_matched.exists():
                        workspace.new_StrmShp_matched.unlink()
                    workspace.DEM_StrmShp = workspace.new_StrmShp_matched
                    return
            
            # Whitebox does not insert the projection into the shapefile, so we need to do that here.
            prj_file = workspace.new_StrmShp.with_suffix('.prj')
            with open(prj_file, 'w') as f:
                f.write(dem_raster.projection)

            if lakes is not None:
                wtbx_gdf = Vector(workspace.new_StrmShp).to_geopandas()
                wtbx_gdf = wtbx_gdf[~wtbx_gdf.geometry.intersects(lakes_gdf.union_all())]  # Remove streams that intersect lakes
                wtbx_gdf.to_file(workspace.new_StrmShp)
            
            wbt.repair_stream_vector_topology(str(workspace.new_StrmShp), str(workspace.new_StrmShp), dist=snap_distance)
            wbt.vector_stream_network_analysis(str(workspace.new_StrmShp), str(workspace.filled_dem), str(workspace.new_StrmShp), snap=snap_distance)

            streams_gdf = _conflate_streams_simple(
                source_gdf=streams_gdf, 
                streams_vector=workspace.new_StrmShp, 
                buffer_distance=buffer_distance
            )
            Vector.save_any_geom(streams_gdf, workspace.new_StrmShp_matched)

            _rasterize_streams(str(workspace.new_stream_raster), str(workspace.fixed_dem), str(workspace.new_StrmShp_matched), attribute=workspace.configs.streamflow_source.upstream_id)

            final_streams = gdal.Open(str(workspace.new_stream_raster)).ReadAsArray()
            channel_mask |= (final_streams > 0)
        else:
            workspace.fldpln_bathymetry = workspace.fixed_dem

    workspace.bathy_water_mask.parent.mkdir(parents=True, exist_ok=True)
    make_channel_mask(channel_mask, str(workspace.fixed_dem), str(workspace.bathy_water_mask))

    stream_info = _create_stream_info_table(
        streams_array=final_streams, 
        streams_gdf=streams_gdf, 
        dem_file=workspace.filled_dem,
        source_id_col=workspace.configs.streamflow_source.upstream_id,
        source_ds_col=workspace.configs.streamflow_source.downstream_id
    )
    stream_info.to_parquet(workspace.stream_info_file, index=False)

    # workspace.assigned_dem = workspace.fixed_dem
    workspace.assigned_dem = workspace.fldpln_bathymetry
    workspace.DEM_StrmShp = workspace.new_StrmShp_matched
    workspace.STRM_File_Clean = workspace.new_stream_raster
    workspace.configs.disable_bathymetry = True

def interpolate_bathymetry(
    bathymetry: np.ndarray,
    channel_mask: np.ndarray,
) -> np.ndarray:
    # 2. Get the coordinates of valid (non-NaN) data and the NaNs
    y, x = np.indices(bathymetry.shape)
    valid_mask = ~np.isnan(bathymetry) & channel_mask
    if not np.any(valid_mask):
        return bathymetry  # No valid points to interpolate from, return the original array

    # Extract coordinates and values of known points
    coords = np.column_stack((y[valid_mask], x[valid_mask]))
    values = bathymetry[valid_mask]


    # Coordinates of the missing values we want to fill
    nan_mask = np.isnan(bathymetry) & channel_mask
    nan_coords = np.column_stack((y[nan_mask], x[nan_mask]))

    # 3. Interpolate using griddata ('linear' or 'cubic')
    interpolated_values = griddata(coords, values, nan_coords, method='linear')

    missing = np.isnan(interpolated_values)
    if np.any(missing):
        interpolated_values[missing] = griddata(
            coords,
            values,
            nan_coords[missing],
            method="nearest",
        )

    # 4. Fill the missing points back into the original raster
    bathymetry[nan_mask] = interpolated_values
    return bathymetry

@njit(cache=True, parallel=True, nogil=True)
def smooth_downhill(
    bathymetry: np.ndarray,
    channel_mask: np.ndarray,
    alpha: float = 2.0,
    iterations: int = 5,
) -> np.ndarray:
    """
    Smooth while giving greater weight to lower neighbors.
    """
    offsets = [
        (-1, -1), (-1, 0), (-1, 1),
        ( 0, -1),          ( 0, 1),
        ( 1, -1), ( 1, 0), ( 1, 1),
    ]
    spatial = np.array([
        [0.7, 1.0, 0.7],
        [1.0, 0.0, 1.0],
        [0.7, 1.0, 0.7],
    ])

    rows, cols = np.where(channel_mask)
    for _ in range(iterations):
        new = bathymetry.copy()

        for i in prange(len(rows)):
            r = rows[i]
            c = cols[i]
            center = bathymetry[r, c]
            if not np.isfinite(center):
                continue

            values = [center]
            weights = [1.0]

            for dr, dc in offsets:
                rr = r + dr
                cc = c + dc

                if (
                    0 <= rr < bathymetry.shape[0]
                    and 0 <= cc < bathymetry.shape[1]
                    and channel_mask[rr, cc]
                    and np.isfinite(bathymetry[rr, cc])
                ):
                    z = bathymetry[rr, cc]

                    w = spatial[dr + 1, dc + 1]
                    if z > center:
                        w *= np.exp(-alpha * (z - center))

                    values.append(z)
                    weights.append(w)

            new[r, c] = np.average(values, weights=weights)

        bathymetry = new

    return bathymetry

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
    channel_mask = burn_streams_into_dem(dem, workspace.DEM_StrmShp, dem_ds, source_gdf, lakes, id_col, ds_col)

    dem = smooth_burned_dem(dem, channel_mask, pbar=False)
    output_ds = gdal.GetDriverByName('GTiff').Create(workspace.fixed_dem, dem_ds.RasterXSize, dem_ds.RasterYSize, 1, gdal.GDT_Float32)
    output_ds.WriteArray(dem)
    output_ds.SetGeoTransform(dem_ds.GetGeoTransform())
    output_ds.SetProjection(dem_ds.GetProjection())
    output_ds = None

    return channel_mask, dem

def build_mask_graph(dem: np.ndarray, mask: np.ndarray) -> nx.Graph:
    G = nx.Graph()

    coords = np.argwhere(mask)
    coord_set = set(map(tuple, coords))

    # Add nodes in bulk
    G.add_nodes_from(
        (
            (r, c),
            {
                "elevation": dem[r, c]
            },
        )
        for r, c in coords
    )

    # Only forward directions (avoid duplicate edges)
    directions = [(0, 1), (1, 0), (1, 1)]
    G.add_edges_from(
        ((r, c), (r + dr, c + dc))
        for r, c in coords
        for dr, dc in directions
        if (r + dr, c + dc) in coord_set
    )
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

def smooth_burned_dem(array: np.ndarray, mask: np.ndarray = None, pbar: bool = True, max_difference: float = 0.5) -> np.ndarray:
    array = array.astype(np.float32, copy=False)
    if mask is None:
        mask = (array % max_difference) == 0

    banks = binary_dilation(mask, structure=np.ones((3, 3), dtype=int)).astype(np.bool_)
    banks[mask] = False  # Banks are the cells adjacent to the mask, but not in the mask

    # Get a distance raster for the mask, where the distance is highest in the center, but smallest at the edges
    distance_to_banks = distance_transform_edt(mask)

    G = build_mask_graph(array, mask)
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
            rows, cols = np.nonzero(
                (array[node[0]-1:node[0]+2, node[1]-1:node[1]+2] <= new_elevation) & 
                banks[node[0]-1:node[0]+2, node[1]-1:node[1]+2]
            )
            if len(rows) > 0:
                array[rows + node[0]-1, cols + node[1]-1] = new_elevation + 0.5

            final_node_elevations[node] = new_elevation

    for node, new_elevation in final_node_elevations.items():
        array[node] = new_elevation

    return array

def burn_streams_into_dem(
        dem: np.ndarray,
        streamlines: np.ndarray,
        dem_ds: gdal.Dataset,
        streams_gdf: gpd.GeoDataFrame,
        lakes: np.ndarray | None,
        id_col: str = 'LINKNO',
        ds_col: str = 'DSLINKNO',
        min_feature_size: int = 5
):
    # Rasterize streamlines on the fly
    mem_ds = gdal.GetDriverByName('MEM').Create('', dem.shape[1], dem.shape[0], 1, gdal.GDT_Int32)
    mem_ds.SetGeoTransform(dem_ds.GetGeoTransform())
    mem_ds.SetProjection(dem_ds.GetProjection())
    stream_ds = ogr.Open(streamlines)
    stream_layer = stream_ds.GetLayer()
    gdal.RasterizeLayer(mem_ds, [1], stream_layer, options=[f"ATTRIBUTE={id_col}"])
    mem_ds.FlushCache()
    streams = mem_ds.ReadAsArray()

    G = nx.from_pandas_edgelist(
        streams_gdf[streams_gdf[ds_col] > 0],
        source=id_col,
        target=ds_col,
        create_using=nx.DiGraph
    )
    channel_mask = (dem % 0.5 == 0)
    # Remove any cells that have less than 5 connected neighbors in the mask, as they are likely to be noise or small artifacts
    structure = np.ones((3, 3), dtype=int)
    labels, _ = label(channel_mask, structure=structure)

    counts = np.bincount(labels.ravel())

    keep = counts > min_feature_size
    keep[0] = False
    channel_mask = keep[labels]

    # Mask out lakes
    if lakes is not None:
        channel_mask &= ~lakes

    # Mask out ocean (where elevation == 0)
    ocean_mask = (dem == 0)
    channel_mask &= ~ocean_mask

    # Mask out dem's no data value
    nodata = dem_ds.GetRasterBand(1).GetNoDataValue()
    if nodata is not None:
        channel_mask &= (dem != nodata)

    labels, num_features = label(channel_mask, structure=structure)
    
    # Buffer the stream raster as a new mask
    channel_border = binary_dilation(streams > 0, structure=np.ones((3, 3), dtype=bool)).astype(np.bool_)
    channel_border &= (streams == 0) & (dem > -9998)

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
                burn_linestring(dem, streams, dem_ds, line, stream_id, G, streams_gdf, channel_mask, labels, local_min)
        elif geom.geom_type == "LineString":
            burn_linestring(dem, streams, dem_ds, geom, stream_id, G, streams_gdf, channel_mask, labels, local_min)
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
        channel_mask: np.ndarray, 
        labels: np.ndarray,
        local_min: np.ndarray):
    if linkno not in G:
        return  # Skip if the linkno is not in the graph
    
    inverse_transform = gdal.InvGeoTransform(dem_ds.GetGeoTransform())

    # Find the first and last points of the linestring in pixel coordinates
    x1, y1 = linestring.coords[0]
    col1, row1 = gdal.ApplyGeoTransform(inverse_transform, x1, y1)
    x2, y2 = linestring.coords[-1]
    col2, row2 = gdal.ApplyGeoTransform(inverse_transform, x2, y2)

    # Determine which end is upstream and which is downstream based on topolgy
    if not (0 <= row1 < dem.shape[0] and 0 <= col1 < dem.shape[1]):
        # Traverse the linestring to find the first point that is within the DEM bounds
        for x, y in linestring.coords:
            col, row = gdal.ApplyGeoTransform(inverse_transform, x, y)
            if 0 <= row < dem.shape[0] and 0 <= col < dem.shape[1]:
                col1, row1 = col, row
                break
        else:
            return  # No valid point found within DEM bounds
        
    if not (0 <= row2 < dem.shape[0] and 0 <= col2 < dem.shape[1]):
        # Traverse the linestring in reverse to find the last point that is within the DEM bounds
        for x, y in reversed(linestring.coords):
            col, row = gdal.ApplyGeoTransform(inverse_transform, x, y)
            if 0 <= row < dem.shape[0] and 0 <= col < dem.shape[1]:
                col2, row2 = col, row
                break
        else:
            return  # No valid point found within DEM bounds
        
    row1, col1 = math.floor(row1), math.floor(col1)
    row2, col2 = math.floor(row2), math.floor(col2)
    row1, col1 = nearest_stream_raster_pixel(streams, row1, col1, linkno)
    row2, col2 = nearest_stream_raster_pixel(streams, row2, col2, linkno)

    # Can't rely on elevation to determine upstream/downstream because the DEM may have been modified by previous burns. (unless no up or downstream nodes exist) Instead, use the topology of the stream network.
    upstream_nodes = list(G.predecessors(linkno))
    downstream_nodes = list(G.successors(linkno))
    if upstream_nodes:          
        # Check if the first point is an endpoint in the upstream node. If not, then the first point is downstream and we need to reverse the line.
        upstream_node = upstream_nodes[0]
        upstream_geom = streams_gdf.at[upstream_node, 'geometry']
        if not Point(x1, y1).touches(upstream_geom):
            linestring = reverse(linestring)
            col1, row1, col2, row2 = col2, row2, col1, row1  # Swap the coordinates as well
    elif downstream_nodes and any(G.successors(downstream_nodes[0])):  # Only check downstream if the downstream node has a successor (i.e., it is not an outlet)
        # Check if the last point is an endpoint in the downstream node. If not, then the last point is upstream and we need to reverse the line.
        downstream_node = downstream_nodes[0]
        downstream_geom = streams_gdf.at[downstream_node, 'geometry']
        if not Point(x2, y2).touches(downstream_geom):
            linestring = reverse(linestring)
            col1, row1, col2, row2 = col2, row2, col1, row1  # Swap the coordinates as well
    else:
        # Fallback to using elevation to determine upstream/downstream if no upstream or downstream nodes exist.
        elev1 = dem[row1, col1]
        elev2 = dem[row2, col2]
        if elev1 < elev2:
            # Reverse the line
            linestring = reverse(linestring)
            col1, row1, col2, row2 = col2, row2, col1, row1  # Swap the coordinates as well

    col2, row2 = col1, row1

    # Check if we have any upstream ids.
    upstream_ids = list(G.predecessors(linkno))
    last_elevation = np.inf
    for upstream_id in upstream_ids:
        us_row, us_col = nearest_stream_raster_pixel(streams, row1, col1, upstream_id)
        if 0 <= us_row < dem.shape[0] and 0 <= us_col < dem.shape[1] and dem[us_row, us_col] < last_elevation:
            last_elevation = dem[us_row, us_col]


    rows = [row1]
    cols = [col1]
    coords = np.asarray(linestring.coords)
    xs = coords[:, 0]
    ys = coords[:, 1]
    for i in range(1, len(linestring.coords)):
        col, row = gdal.ApplyGeoTransform(inverse_transform, xs[i], ys[i])
        row, col = math.floor(row), math.floor(col)
        row, col = nearest_stream_raster_pixel(streams, row, col, linkno)
        if not (0 <= row < dem.shape[0] and 0 <= col < dem.shape[1]):
            continue

        if row != rows[-1] or col != cols[-1]:
            rows.append(row)
            cols.append(col)

    started_in_mask = channel_mask[row1, col1]
    write = True
    will_come_back = True

    for i, (row1, col1, row2, col2) in enumerate(zip(rows[:-1], cols[:-1], rows[1:], cols[1:])):
        if started_in_mask and write and not channel_mask[row2, col2]:
            if will_come_back:
                # Check if this linestring will eventually go back in the mask (and that point is connected to this one)
                for j in range(i + 1, len(rows)):
                    row3, col3 = rows[j], cols[j]
                    if channel_mask[row3, col3] and labels[row1, col1] == labels[row3, col3]:
                        # We can stop burning because we know we will come back (the source stream wandered away slightly from DEM channel)
                        write = False
                        break
                else:
                    will_come_back = False
        elif started_in_mask and not write and channel_mask[row2, col2]:
            write = True
        elif not started_in_mask and channel_mask[row2, col2]:
            started_in_mask = True

        if not write:
            continue

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


@dataclass
class NodeSig:
    kind: NodeType                      # outlet, confluence, source, interior
    local_order: int


@dataclass
class ReachSig:
    geom: LineString
    buffer_geom: Polygon
    centroid: Point


def build_graph(gdf: gpd.GeoDataFrame, id_col: str, ds_col: str) -> nx.DiGraph:
    """
    Build a downstream-directed DAG where each graph node is a stream reach.
    Node attributes store the reach geometry and derived properties.
    The downstream reach is represented by a directed edge: upstream_reach -> downstream_reach.
    """
    G = nx.DiGraph()

    if gdf.empty:
        return G

    if "geometry" not in gdf.columns:
        raise ValueError("GeoDataFrame must have a geometry column.")

    for row in gdf.itertuples(index=False):
        rid = getattr(row, id_col)
        geom = getattr(row, "geometry")
        buffered = getattr(row, "buffered")
        if geom is None or geom.is_empty:
            continue

        if isinstance(geom, GeometryCollection):
            geom = max((g for g in geom.geoms if isinstance(g, LineString)), key=lambda x: x.length, default=None)
            if geom is None:
                continue

        multi_geom = None

        if not isinstance(geom, LineString):
            # if a MultiLineString slipped through, merge if possible
            if isinstance(geom, MultiLineString):
                multi_geom = geom
                merged = line_merge(geom)
                if isinstance(merged, LineString):
                    geom = merged
                elif isinstance(merged, MultiLineString):
                    # The downstream is the geom which tocuhes the gdf, without this id in it
                    touching = gdf[(gdf[id_col] != rid) & gdf.touches(merged)]
                    if touching.empty:
                        continue
                    downstream = touching[
                        ~touching[ds_col].isin(touching[id_col])
                    ]
                    geoms = [g for g in merged.geoms if g.geom_type == "LineString"]
                    downstream_geom = [g for g in geoms if any(g.touches(d) for d in downstream.geometry)]
                    if downstream_geom:
                        geom = downstream_geom[0] # We want this geometry to match with wtbx
                    else:
                        # Just use the longest geom if no downstream geom is found
                        geom = max(geoms, key=lambda x: x.length)
                else:
                    continue
            else:
                continue

        ds_id = getattr(row, ds_col, None)
        if pd_is_na(ds_id):
            ds_id = 0

        G.add_node(
            rid,
            geometry=geom,
            multi_geom=multi_geom,
            centroid=geom.centroid,
            STRAHLER=getattr(row, "STRAHLER", None),
            buffered=buffered
        )

    for row in gdf.itertuples(index=False):
        rid = getattr(row, id_col)
        ds_id = getattr(row, ds_col, None)
        if pd_is_na(ds_id):
            continue
        try:
            ds_id = int(ds_id)
        except Exception:
            # allow string ids too
            pass
        if ds_id in (0, -1, None):
            continue
        if rid in G and ds_id in G:
            G.add_edge(rid, ds_id)

    return G


# ---------------------------------------------------------------------
# Signatures
# ---------------------------------------------------------------------

def compute_node_signatures(G: nx.DiGraph) -> dict[str, NodeSig]:
    sig: dict[str, NodeSig] = {}
    for n, data in G.nodes(data=True):
        sig[n] = NodeSig(
            kind=classify_node(G, n, data['STRAHLER']),
            local_order=data['STRAHLER'],
        )
    return sig


def compute_reach_signatures(G: nx.DiGraph) -> dict[str, ReachSig]:
    """
    In this graph representation, every graph node is a reach.
    So the reach signature is attached to each node id.
    """
    sig: dict[str, ReachSig] = {}
    for n, data in G.nodes(data=True):
        geom = data["geometry"]
        prepare(geom)
        sig[n] = ReachSig(
            geom=geom,
            buffer_geom=data['buffered'],
            centroid=data['centroid'],
        )
    return sig

class NodeType(Enum):
    ISOLATED = 1
    OUTLET = 2
    SOURCE = 3
    MERGED_SOURCE = 4
    CONFLUENCE = 5
    INTERIOR = 6

    def is_outlet(self) -> bool:
        return self == NodeType.OUTLET or self == NodeType.ISOLATED
    
    def is_source(self) -> bool:
        return self == NodeType.SOURCE or self == NodeType.MERGED_SOURCE or self == NodeType.ISOLATED

def classify_node(G: nx.DiGraph, n: str, local_order: int) -> NodeType:
    indeg = G.in_degree(n)
    outdeg = G.out_degree(n)
    if indeg == 0 and outdeg == 0:
        return NodeType.ISOLATED
    if outdeg == 0:
        return NodeType.OUTLET
    if indeg == 0:
        if local_order > 1:
            return NodeType.MERGED_SOURCE
        return NodeType.SOURCE
    if indeg > 1:
        return NodeType.CONFLUENCE
    return NodeType.INTERIOR


def strm_order_similarity(A: NodeSig, B: NodeSig) -> float:
    if A.local_order is None or B.local_order is None:
        return 0.0
    d = abs(A.local_order - B.local_order)
    return math.exp(-d/5)

def point_similarity(p1: Point, p2: Point, scale: float) -> float:
    d = p1.distance(p2)
    return math.exp(-d / max(scale, 1.0))

def pd_is_na(x) -> bool:
    try:
        return bool(x is None or (isinstance(x, float) and np.isnan(x)))
    except Exception:
        return x is None

# ---------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------

def downstream_centroid(G: nx.DiGraph, node: int, reach_sig: dict[str, ReachSig]) -> Point:
    if G.out_degree(node) > 0:
        downstream = list(G.successors(node))[0]
        # Find the endpoint which touches the downstream node
        geom = reach_sig[node].geom
        start_pt = Point(geom.coords[0])
        end_pt = Point(geom.coords[-1])
        if start_pt.touches(reach_sig[downstream].geom):
            return start_pt
        return end_pt

    if G.in_degree(node) > 0:
        upstream = list(G.predecessors(node))[0]
        geom = reach_sig[node].geom
        start_pt = Point(geom.coords[0])
        end_pt = Point(geom.coords[-1])
        if start_pt.touches(reach_sig[upstream].geom):
            return end_pt
        return start_pt
    
    # Just return the centroid if isolated
    return reach_sig[node].centroid

def upstream_centroid(G: nx.DiGraph, node: int, reach_sig: dict[str, ReachSig]) -> Point:
    if G.in_degree(node) > 0:
        upstream = list(G.predecessors(node))[0]
        # Find the endpoint which touches the upstream node
        geom = reach_sig[node].geom
        start_pt = Point(geom.coords[0])
        end_pt = Point(geom.coords[-1])
        if start_pt.touches(reach_sig[upstream].geom):
            return start_pt
        return end_pt

    if G.out_degree(node) > 0:
        downstream = list(G.successors(node))[0]
        geom = reach_sig[node].geom
        start_pt = Point(geom.coords[0])
        end_pt = Point(geom.coords[-1])
        if end_pt.touches(reach_sig[downstream].geom):
            return start_pt
        return end_pt
    
    # Assume first point is upstream if isolated
    return Point(reach_sig[node].geom.coords[0])

def _conflate_streams_simple(
        source_gdf: gpd.GeoDataFrame,
        streams_vector: os.PathLike,
        buffer_distance: float,
) -> gpd.GeoDataFrame:
    """
    We assume that source and wtbx streams are super closely aligned already.
    We can just iterate over each source stream, find the one that intersects it the most, and assign that
    """
    wtbx_vector = Vector(streams_vector)
    wtbx_gdf = wtbx_vector.to_geopandas()

    source_gdf['corresponding_wtbx_id'] = None
    wtbx_sindex = wtbx_gdf.sindex

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, message="Geometry is in a geographic CRS")
        wtbx_gdf['buffered'] = wtbx_gdf.buffer(buffer_distance, resolution=2)
        source_gdf['buffered'] = source_gdf.buffer(buffer_distance, resolution=2)

    for source_row in source_gdf.itertuples():
        source_geom: Polygon = source_row.buffered
        possible_matches_index = list(wtbx_sindex.intersection(source_geom.bounds))
        possible_matches: gpd.GeoDataFrame = wtbx_gdf.iloc[possible_matches_index]
        best_match = None
        best_area = 0
        for wtbx_row in possible_matches.itertuples():
            intersection_area = source_geom.intersection(wtbx_row.buffered).area
            if intersection_area > best_area:
                best_area = intersection_area
                best_match = wtbx_row

        if best_match is not None:
            source_gdf.at[source_row.Index, 'corresponding_wtbx_id'] = best_match.Index

    # Transfer all the attributes from the source to the wtbx gdf
    source_gdf["corresponding_wtbx_id"] = source_gdf["corresponding_wtbx_id"].astype('Int64')
    wtbx_gdf = wtbx_gdf[['FID']].merge(
        source_gdf.drop(columns=['buffered', 'FID']), 
        left_index=True, 
        right_on='corresponding_wtbx_id'
    ).drop(columns='corresponding_wtbx_id')

    wtbx_gdf = gpd.GeoDataFrame(wtbx_gdf, geometry='geometry', crs=source_gdf.crs)

    return wtbx_gdf


def is_linear_reach(G: nx.DiGraph, upstream: int, downstream: int):
    path = nx.shortest_path(G, upstream, downstream)

    for node in path[1:]:
        if G.in_degree(node) != 1:
            return False

    return True

def _conflate_streams(
    source_gdf: gpd.GeoDataFrame,
    streams_vector: os.PathLike,
    buffer_distance: float,
    dem_proj: str,
    dem_bbox: tuple[float, float, float, float],
    source_id_col: str = "LINKNO",
    source_ds_col: str = "DSLINKNO",
    wtbx_id_col: str = "FID",
    wtbx_ds_col: str = "DS_LINK_ID",
    strm_order_col: str = "strmOrder",
) -> gpd.GeoDataFrame:
    """
    Read files, build graphs, match the stream networks, and return the mapping state.
    """
    wtbx_vector = Vector(streams_vector)
    wtbx_gdf = wtbx_vector.to_geopandas().set_crs(dem_proj)
    bounds = wtbx_vector.epsg_4326_bbox

    # Filter source_gdf to only include streams that intersect the bounds of the wtbx_gdf, shrunk by a small distance
    source_gdf = source_gdf[source_gdf.intersects(box(*bounds))].copy()
    source_gdf = source_gdf.to_crs(dem_proj)
    min_strm_order = source_gdf[strm_order_col].min()

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, message="Geometry is in a geographic CRS")
        wtbx_gdf['buffered'] = wtbx_gdf.buffer(buffer_distance, resolution=2)
        source_gdf['buffered'] = source_gdf.simplify(tolerance=1e-5).buffer(buffer_distance, resolution=2)

    GA = build_graph(source_gdf, id_col=source_id_col, ds_col=source_ds_col)
    GB = build_graph(wtbx_gdf, id_col=wtbx_id_col, ds_col=wtbx_ds_col)

    assert nx.is_directed_acyclic_graph(GA), "GA must be a DAG."
    assert nx.is_directed_acyclic_graph(GB), "GB must be a DAG."

    for node in nx.topological_sort(GA):
        orders = [
            GA.nodes[p]['STRAHLER']
            for p in GA.predecessors(node)
        ]

        if len(orders) == 0:
            GA.nodes[node]['STRAHLER'] = min_strm_order # The source gdf has headwater streams merged, so there is no stream order 1. I assume that we will always have a headwater stream in the AOI for this to work.
            continue

        max_order = max(orders)

        if orders.count(max_order) >= 2:
            GA.nodes[node]['STRAHLER'] = max_order + 1
        else:
            GA.nodes[node]['STRAHLER'] = max_order

    for i, node in enumerate(nx.topological_sort(GB)):
        GB.nodes[node]['topological_order'] = i

    # Use a projected CRS for distance/length calculations
    wtbx_gdf[source_id_col] = np.nan

    A_reach_sig = compute_reach_signatures(GA)
    B_reach_sig = compute_reach_signatures(GB)

    A_node_sig = compute_node_signatures(GA)
    B_node_sig = compute_node_signatures(GB)

    # Get outlets of source graph
    source_headwaters = {n for n, sig in A_node_sig.items() if sig.kind.is_source()}
    wtbx_headwaters = {n for n, sig in B_node_sig.items() if sig.kind.is_source()}
    wtbx_potential_headwaters = {
        n for n, sig in B_node_sig.items() if sig.local_order == min_strm_order and all(B_node_sig[pred].local_order < min_strm_order for pred in GB.predecessors(n)) # The predecessors are both headwaters
    }
    wtbx_potential_headwaters -= wtbx_headwaters
    non_source_headwaters = set(GA.nodes()) - source_headwaters

    # Match headwaters
    best_source_scores = {}
    scores = {}
    for source_headwater in source_headwaters:
        best_score = 0
        best_wtbx = None
        source_centroid = A_reach_sig[source_headwater].centroid
        source_buffer = A_reach_sig[source_headwater].buffer_geom

        for wtbx_outlet in sorted(wtbx_headwaters, key=lambda x: B_reach_sig[x].geom.distance(source_centroid))[:5]:
            wtbx_buffer = B_reach_sig[wtbx_outlet].buffer_geom
            intersection_area = source_buffer.intersection(wtbx_buffer).area
            score = intersection_area

            if score > best_score:
                best_score = score
                best_wtbx = wtbx_outlet

        for wtbx_outlet in sorted(wtbx_potential_headwaters, key=lambda x: B_reach_sig[x].geom.distance(source_centroid))[:5]:
            wtbx_buffer = B_reach_sig[wtbx_outlet].buffer_geom
            intersection_area = source_buffer.intersection(wtbx_buffer).area
            score = intersection_area
            if score > best_score:
                best_score = score
                best_wtbx = wtbx_outlet

        if best_wtbx is not None and best_score > scores.get(best_wtbx, 0):
            best_source_scores[source_headwater] = (best_wtbx, source_headwater)
            scores[best_wtbx] = best_score

    # Check the computed wtbx headwaters, and see if there is a better match for any of them. If so, replace the match.
    for source_headwater, (wtbx_headwater, _) in list(best_source_scores.items()):
        best_score = scores.get(wtbx_headwater, 0)
        best_source = None
        wtbx_centroid = B_reach_sig[wtbx_headwater].centroid
        wtbx_buffer = B_reach_sig[wtbx_headwater].buffer_geom

        for source_headwater2 in sorted(non_source_headwaters, key=lambda x: A_reach_sig[x].geom.distance(wtbx_centroid))[:5]:
            source_buffer2 = A_reach_sig[source_headwater2].buffer_geom
            intersection_area = source_buffer2.intersection(wtbx_buffer).area
            score = intersection_area

            if score > best_score:
                best_score = score
                best_source = source_headwater2

        if best_source is not None and best_source != source_headwater:
            if not nx.has_path(GA, best_source, source_headwater) and nx.has_path(GA, source_headwater, best_source) and is_linear_reach(GA, source_headwater, best_source):
                # We matched better with a downstream segment, but let us actually keep the upstream one instead (in a linear reach, i.e., no upstream tributaries)
                continue
            del best_source_scores[source_headwater]
            best_source_scores[best_source] = (wtbx_headwater, best_source)

    headwater_fids = set(f for f, _ in best_source_scores.values())
    headwater_fid_to_linkno = {f: l for f, l in best_source_scores.values()}

    headwaters_with_no_outlet = headwater_fids
    headwater_fids_accounted_for = set()
    best_outlet_linkno_to_fid = {}
    if headwaters_with_no_outlet:
        visited = set()
        # For each headwater, find the corresponding outlet
        for headwater_fid in headwaters_with_no_outlet:
            if headwater_fid in visited:
                continue

            headwater_linkno = headwater_fid_to_linkno[headwater_fid]
            # Find the outlet for this headwater
            outlet_fid = None
            descendants = list(nx.descendants(GB, headwater_fid))
            if len(descendants) == 0:
                # This headwater is an outlet, so we can just use it as the outlet
                headwater_fids_accounted_for.add(headwater_fid)
                best_outlet_linkno_to_fid[headwater_linkno] = headwater_fid
                continue
            else:
                for fid in nx.descendants(GB, headwater_fid):
                    if GB.out_degree(fid) == 0:
                        outlet_fid = fid
                        break

            if outlet_fid is None:
                # Remove these headwaters from the best_source_scores, since they don't have a corresponding outlet
                del best_source_scores[headwater_fid_to_linkno[headwater_fid]]
                visited.add(headwater_fid)
                continue

            # Find all headwater fids that connect to this outlet, with no current outlet
            outlet_headwaters = set(f for f in nx.ancestors(GB, outlet_fid)) & headwaters_with_no_outlet
            
            # Find the most common downstream linkno for these headwaters in GA
            headwater_linknos = [headwater_fid_to_linkno[f] for f in outlet_headwaters]
            current_linkno_outlet = headwater_linknos[0]
            found = False
            while True:
                try:
                    current_linkno_outlet = next(GA.successors(current_linkno_outlet))
                except StopIteration:
                    break
                if all(nx.has_path(GA, hdwtr, current_linkno_outlet) for hdwtr in headwater_linknos):
                    found = True
                    break

            if not found:
                # Remove these headwaters from the best_source_scores, since they don't have a corresponding outlet
                for headwater_fid in outlet_headwaters:
                    del best_source_scores[headwater_fid_to_linkno[headwater_fid]]
                    visited.add(headwater_fid)
                continue

            headwater_fids_accounted_for.update(outlet_headwaters)
            best_outlet_linkno_to_fid[current_linkno_outlet] = outlet_fid

        headwater_fids = headwater_fids_accounted_for

    pairs = list((fid, linkno) for linkno, fid in best_outlet_linkno_to_fid.items()) + list(best_source_scores.values())
    final_matches = {}
    linkno_to_fid = defaultdict(set)
    for fid, linkno in pairs:
        final_matches[fid] = linkno
        linkno_to_fid[linkno].add(fid)

    for outlet_linkno, outlet_fid in best_outlet_linkno_to_fid.items():
        allowable_linknos = defaultdict(set)
        watershed_fids = set(nx.ancestors(GB, outlet_fid)) | {outlet_fid}
        # Find common downstream segment for each headwater
        headwaters_of_outlet_that_are_assigned = watershed_fids & headwater_fids
        for fid_source_1, fid_source_2 in combinations(headwaters_of_outlet_that_are_assigned, 2):
            if fid_source_1 == fid_source_2:
                continue

            # Find the first fid which has both as upstream.
            # Do the same for linkno
            confluence_fid = fid_source_1
            while True:
                if nx.has_path(GB, fid_source_2, confluence_fid):
                    break
                confluence_fid = next(GB.successors(confluence_fid))

            confluence_linkno = headwater_fid_to_linkno[fid_source_1]
            do_these_headwaters_connect = True
            while True:
                if nx.has_path(GA, headwater_fid_to_linkno[fid_source_2], confluence_linkno):
                    break
                try:
                    confluence_linkno = next(GA.successors(confluence_linkno))
                except StopIteration:
                    do_these_headwaters_connect = False
                    break
            if not do_these_headwaters_connect:
                continue

            if confluence_fid in final_matches and final_matches[confluence_fid] != confluence_linkno:
                # Choose the most downstream linkno
                if nx.has_path(GA, confluence_linkno, final_matches[confluence_fid]):
                    confluence_linkno = final_matches[confluence_fid]

            if confluence_linkno in linkno_to_fid:
                # Fill all the fids between the linkno and the outlet with the same linkno
                other_fids = linkno_to_fid[confluence_linkno]
                to_add = set()
                for fid in other_fids:
                    if nx.has_path(GB, confluence_fid, fid):
                        for fid in nx.shortest_path(GB, confluence_fid, fid):
                            final_matches[fid] = confluence_linkno
                            to_add.add(fid)
                    elif nx.has_path(GB, fid, confluence_fid):
                        for fid in nx.shortest_path(GB, fid, confluence_fid):
                            final_matches[fid] = confluence_linkno
                            to_add.add(fid)
                    else:
                        continue

                linkno_to_fid[confluence_linkno].update(to_add)

            final_matches[confluence_fid] = confluence_linkno
            linkno_to_fid[confluence_linkno].add(confluence_fid)

            source1_path = set(nx.shortest_path(GA, headwater_fid_to_linkno[fid_source_1], confluence_linkno))
            for fid in nx.shortest_path(GB, fid_source_1, confluence_fid):
                allowable_linknos[fid].update(source1_path)

            source2_path = set(nx.shortest_path(GA, headwater_fid_to_linkno[fid_source_2], confluence_linkno))
            for fid in nx.shortest_path(GB, fid_source_2, confluence_fid):
                allowable_linknos[fid].update(source2_path)

            try:
                source3_path = set(nx.shortest_path(GA, confluence_linkno, outlet_linkno))
            except nx.NetworkXNoPath:
                source3_path = set() # Occurs if outlet_fid == confluence_fid
            for fid in nx.shortest_path(GB, confluence_fid, outlet_fid):
                allowable_linknos[fid].update(source3_path)

        if len(headwaters_of_outlet_that_are_assigned) == 1:
            source_path = set(nx.shortest_path(GA, headwater_fid_to_linkno[next(iter(headwaters_of_outlet_that_are_assigned))], outlet_linkno))
            for fid in nx.shortest_path(GB, next(iter(headwaters_of_outlet_that_are_assigned)), outlet_fid):
                allowable_linknos[fid].update(source_path)

        for fid, allowable in allowable_linknos.items():
            if fid in final_matches:
                continue

            new_allowable = []
            for linkno in allowable:
                existing_fids = linkno_to_fid[linkno]
                # Only add if the existing fids are immediately up or downstream of this fid
                if existing_fids:
                    S = GB.subgraph(existing_fids | {fid})
                    if not nx.is_weakly_connected(S):
                        continue

                    # Ensure in and out degree all <= 1
                    if any(S.in_degree(n) > 1 or S.out_degree(n) > 1 for n in S.nodes):
                        continue
                    
                new_allowable.append(linkno)

            allowable = new_allowable
            allowable = sorted(allowable, key=lambda x: A_reach_sig[x].geom.distance(B_reach_sig[fid].geom))[:4]

            fid_geom = B_reach_sig[fid].geom
            fid_buffer = B_reach_sig[fid].buffer_geom
            best_score = 0
            best_linkno = None
            for linkno in allowable:
                linkno_buffer = A_reach_sig[linkno].buffer_geom
                score = fid_buffer.intersection(linkno_buffer).area
                if score > best_score:
                    best_score = score
                    best_linkno = linkno

            if best_linkno is None:
                # Use hausdorff distance to find the closest linkno
                fid_geom = B_reach_sig[fid].geom
                best_score = np.inf
                for linkno in allowable:
                    linkno_geom = A_reach_sig[linkno].geom
                    score = fid_geom.hausdorff_distance(linkno_geom)
                    if score < best_score:
                        best_score = score
                        best_linkno = linkno

            final_matches[fid] = best_linkno
            linkno_to_fid[best_linkno].add(fid)


    # Finally assign to the gdf
    wtbx_gdf[source_id_col] = wtbx_gdf['FID'].map(lambda x: final_matches.get(x, np.nan))
    wtbx_gdf = wtbx_gdf[wtbx_gdf[source_id_col].notna()].copy()

    min_bounds = (
        dem_bbox[0] + buffer_distance,
        dem_bbox[1] + buffer_distance,
        dem_bbox[2] - buffer_distance,
        dem_bbox[3] - buffer_distance
    )
    geom_intersects_bounds = wtbx_gdf.geometry.intersects(box(*min_bounds).boundary)
    if not geom_intersects_bounds.all():
        wtbx_gdf = wtbx_gdf[~wtbx_gdf.geometry.intersects(box(*min_bounds).boundary)].copy() # Remove streams that touch near the border of the AOI, since they are likely to be the most incorrect
    wtbx_gdf['topological_order'] = wtbx_gdf['FID'].map(lambda x: GB.nodes[x]['topological_order'])

    wtbx_gdf = wtbx_gdf.dissolve(source_id_col, as_index=False, sort=False, aggfunc={
        'FID': partial(select_most_downstream_fid, GB=GB),
        'topological_order': 'min',
        'STRAHLER': 'max',
    })
    wtbx_gdf['geometry'] = wtbx_gdf['geometry'].apply(lambda geom: line_merge(geom))
    wtbx_gdf[source_id_col] = wtbx_gdf[source_id_col].astype(int)

    # Compute DSLINKNO for the conflated streams
    fid_to_dsfid = {}
    for fid in wtbx_gdf['FID']:
        downstream_fid = next(GB.successors(fid), None)
        if downstream_fid is not None:
            fid_to_dsfid[fid] = downstream_fid
        else:
            fid_to_dsfid[fid] = -1

    wtbx_gdf[source_ds_col] = wtbx_gdf['FID'].map(fid_to_dsfid)
    wtbx_gdf[source_ds_col] = wtbx_gdf[source_ds_col].map(lambda x: final_matches.get(x, -1) if x != -1 else -1)
    wtbx_gdf.loc[wtbx_gdf[source_id_col] == wtbx_gdf[source_ds_col], source_ds_col] = -1 # Remove self-loops (since we filtered out the borders)

    # Move linkno to the front
    cols = wtbx_gdf.columns.tolist()
    cols.insert(0, cols.pop(cols.index(source_id_col)))
    cols.insert(1, cols.pop(cols.index(source_ds_col)))
    wtbx_gdf = wtbx_gdf[cols]
    wtbx_gdf = wtbx_gdf[wtbx_gdf['geometry'].geom_type == 'LineString'].copy()

    return wtbx_gdf

def select_most_downstream_fid(gdf: gpd.GeoSeries, GB: nx.DiGraph) -> gpd.GeoDataFrame:
    S = GB.subgraph(gdf.values)
    return [s for s in S.nodes if S.out_degree(s) == 0][0]


def _rasterize_streams(stream_raster: str, dem: str, streams_vector: str, attribute: str):
    dem_ds: gdal.Dataset = gdal.Open(dem)
    stream_ds: gdal.Dataset = gdal.GetDriverByName("GTiff").Create(stream_raster, dem_ds.RasterXSize, dem_ds.RasterYSize, 1, gdal.GDT_Int32, ['COMPRESS=DEFLATE', 'PREDICTOR=2'])
    stream_ds.SetGeoTransform(dem_ds.GetGeoTransform())
    stream_ds.SetProjection(dem_ds.GetProjection())
    options = gdal.RasterizeOptions(attribute=attribute)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        gdal.Rasterize(stream_ds, streams_vector, options=options)

def _create_stream_info_table(
        streams_array: np.ndarray, 
        streams_gdf: gpd.GeoDataFrame,
        dem_file: str,
        source_id_col: str,
        source_ds_col: str
        ) -> pd.DataFrame:
    ds: gdal.Dataset = gdal.Open(dem_file)
    dem: np.ndarray = ds.ReadAsArray()
    gt = ds.GetGeoTransform()
    nrows, ncols = streams_array.shape
    
    G: nx.DiGraph = nx.from_pandas_edgelist(
        streams_gdf[streams_gdf[source_ds_col] > 0],
        source=source_id_col,
        target=source_ds_col,
        create_using=nx.DiGraph()
    )

    values, counts = np.unique(streams_array, return_counts=True)
    linkno_counts = dict(zip(values, counts))

    output_table = []
    for idx, row in streams_gdf.iterrows():
        line: LineString = row.geometry
        linkno = row[source_id_col]

        # Get the pixel coordinates of the start and end points
        start_point = Point(line.coords[0])
        end_point = Point(line.coords[-1])
        start_col_raw = (start_point.x - gt[0]) / gt[1]
        start_row_raw = (start_point.y - gt[3]) / gt[5]
        start_col = round(start_col_raw)
        start_row = round(start_row_raw)
        end_col_raw = (end_point.x - gt[0]) / gt[1]
        end_row_raw = (end_point.y - gt[3]) / gt[5]
        end_col = round(end_col_raw)
        end_row = round(end_row_raw)

        # For both the start and end points, check if we are right on the stream pixel. If not, check the neighbors and choose whichever is closest to the original point. This is to account for slight misalignments between the stream vector and raster.
        if not (0 <= start_row < streams_array.shape[0] and 0 <= start_col < streams_array.shape[1]) or streams_array[start_row, start_col] != linkno:
            min_dist = float('inf')
            closest_row, closest_col = None, None
            for dr in range(-1, 2):
                for dc in range(-1, 2):
                    r = start_row + dr
                    c = start_col + dc
                    if 0 <= r < streams_array.shape[0] and 0 <= c < streams_array.shape[1]:
                        if streams_array[r, c] == linkno:
                            dist = np.sqrt((c - start_col_raw) ** 2 + (r - start_row_raw) ** 2)
                            if dist < min_dist:
                                min_dist = dist
                                closest_row, closest_col = r, c
            if closest_row is not None and closest_col is not None:
                start_row, start_col = closest_row, closest_col

        if not (0 <= end_row < streams_array.shape[0] and 0 <= end_col < streams_array.shape[1]) or streams_array[end_row, end_col] != linkno:
            min_dist = float('inf')
            closest_row, closest_col = None, None
            for dr in range(-1, 2):
                for dc in range(-1, 2):
                    r = end_row + dr
                    c = end_col + dc
                    if 0 <= r < streams_array.shape[0] and 0 <= c < streams_array.shape[1]:
                        if streams_array[r, c] == linkno:
                            dist = np.sqrt((c - end_col_raw) ** 2 + (r - end_row_raw) ** 2)
                            if dist < min_dist:
                                min_dist = dist
                                closest_row, closest_col = r, c
            if closest_row is not None and closest_col is not None:
                end_row, end_col = closest_row, closest_col

        if linkno in G:
            upstream_nodes = list(G.predecessors(linkno))
            downstream_nodes = list(G.successors(linkno))
        else:
            upstream_nodes = []
            downstream_nodes = []
        if upstream_nodes:          
            # Check if the first point is an endpoint in the upstream node. If not, then the first point is downstream and we need to reverse the line.
            upstream_node = upstream_nodes[0]
            upstream_geom = streams_gdf.loc[streams_gdf[source_id_col] == upstream_node, 'geometry'].values[0]
            if not start_point.touches(upstream_geom):
                start_row, start_col, end_row, end_col = end_row, end_col, start_row, start_col
        elif downstream_nodes and any(G.successors(downstream_nodes[0])):  # Only check downstream if the downstream node has a successor (i.e., it is not an outlet)
            # Check if the last point is an endpoint in the downstream node. If not, then the last point is upstream and we need to reverse the line.
            downstream_node = downstream_nodes[0]
            downstream_geom = streams_gdf.loc[streams_gdf[source_id_col] == downstream_node, 'geometry'].values[0]
            if not end_point.touches(downstream_geom):
                start_row, start_col, end_row, end_col = end_row, end_col, start_row, start_col
        else:
            # Fallback to using elevation to determine upstream/downstream if no upstream or downstream nodes exist.
            elev1 = dem[start_row, start_col]
            elev2 = dem[end_row, end_col]
            if elev1 < elev2:
                # Reverse the line
                start_row, start_col, end_row, end_col = end_row, end_col, start_row, start_col

        # Length is simply the number of pixels with that linkno
        length = linkno_counts.get(linkno, 0)

        # convert start_row, start_col, end_row, end_col to single value
        start_idx = start_row * ncols + start_col
        end_idx = end_row * ncols + end_col

        output_table.append((start_idx, end_idx, length, linkno))

    stream_info = pd.DataFrame(output_table, 
                 columns=['start_pixel', 'end_pixel', 'length', 'stream_id'])
    return stream_info

def make_channel_mask(channel_mask: np.ndarray, dem: str, water_mask: str):
    dem_ds: gdal.Dataset = gdal.Open(dem)
    out_ds: gdal.Dataset = gdal.GetDriverByName('GTiff').Create(water_mask, dem_ds.RasterXSize, dem_ds.RasterYSize, 1, gdal.GDT_Byte)
    out_ds.WriteArray(channel_mask)
    out_ds.SetGeoTransform(dem_ds.GetGeoTransform())
    out_ds.SetProjection(dem_ds.GetProjection())
    out_ds = None