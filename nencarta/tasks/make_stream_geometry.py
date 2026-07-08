import json
from pathlib import Path

import pandas as pd
import geopandas as gpd

from nencarta.logger import LOG
from nencarta.api.raster import Raster
from nencarta.api.vector import Vector
from nencarta.workspace import Workspace
from nencarta.api.configs import NencartaConfig

def _filter_streams_with_lake_json(lake_filter_json: str, stream_df: gpd.GeoDataFrame, rivid_field: str) -> gpd.GeoDataFrame:
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

def _filter_streams_by_stream_order(stream_df: gpd.GeoDataFrame, strm_order_field: str, strm_order_low: float | None, strm_order_high: float | None) -> gpd.GeoDataFrame:
    if strm_order_field not in stream_df.columns:
        LOG.warning(f"StrmOrder_Field '{strm_order_field}' not found in stream shapefile; skipping stream order filter.")
        return stream_df

    stream_df[strm_order_field] = pd.to_numeric(stream_df[strm_order_field], errors="coerce")
    if strm_order_low is None:
        strm_order_low = stream_df[strm_order_field].min()
    if strm_order_high is None:
        strm_order_high = stream_df[strm_order_field].max()

    return stream_df[stream_df[strm_order_field].between(strm_order_low, strm_order_high)]


def make_stream_geometry(workspace: Workspace,) -> Path:
    configs = workspace.configs
    if workspace.DEM_StrmShp.exists() and not configs.process_stream_network:
        LOG.info(f"{workspace.DEM_StrmShp} already exists and we aren't making it again...")
        return workspace.DEM_StrmShp
    
    bbox = Raster(workspace.assigned_dem).epsg_4326_bbox

    if configs.flowline:
        LOG.info(f"Using provided flowline {configs.flowline} for stream geometry.")
        streamlines = [configs.flowline]
    elif configs.source_flowlines:
        LOG.info(f"Getting streamlines in extent {bbox}...")
        streamlines = Vector.get_streamlines_in_extent(bbox, configs.source_flowlines)
        if not streamlines:
            if configs.raise_errors_if_nothing_in_domain:
                raise ValueError("No RFS stream geometry files intersect the DEM extent.")
            
            return None

    if len(streamlines) == 1:
        gdf = Vector(streamlines[0]).to_geopandas(bbox_epsg_4326=bbox)
    else:
        gdf = pd.concat([Vector(path).to_geopandas(bbox_epsg_4326=bbox) for path in streamlines], ignore_index=True, copy=False)

    gdf = gdf[~gdf.geometry.isna()]

    if "LINKNO" not in gdf.columns and "COMID" in gdf.columns:
        gdf["LINKNO"] = gdf["COMID"]
    if "COMID" not in gdf.columns and "LINKNO" in gdf.columns:
        gdf["COMID"] = gdf["LINKNO"]

    if configs.lake_filter_json:
        gdf = _filter_streams_with_lake_json(configs.lake_filter_json, gdf, configs.stream_id_field)

    # if the StrmOrder_Field and StrmOrder_Lower or StrmOrder_Upper are not None use these to filter the StrmShp_gdf
    if configs.StrmOrder_Field and (configs.StrmOrder_Lower is not None or configs.StrmOrder_Upper is not None):
        gdf = _filter_streams_by_stream_order(gdf, configs.StrmOrder_Field, configs.StrmOrder_Lower, configs.StrmOrder_Upper)

    # TODO
    if configs.exclude:
        gdf = gdf[~gdf[configs.stream_id_field].isin(configs.exclude)]

    if gdf.empty:
        if configs.raise_errors_if_nothing_in_domain:
            raise ValueError("No stream geometries intersect the DEM extent.")
        
        return None
    
    kwargs = {}
    if workspace.DEM_StrmShp.suffix.lower().endswith('.parquet'):
        kwargs['compression'] = 'brotli'
        kwargs['write_covering_bbox'] = True

    workspace.DEM_StrmShp.parent.mkdir(parents=True, exist_ok=True)

    Vector.save_any_geom(gdf, workspace.DEM_StrmShp, **kwargs)

    return workspace.DEM_StrmShp