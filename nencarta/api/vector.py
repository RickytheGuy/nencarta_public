from __future__ import annotations

import os
import json
import warnings
from typing import Any


import pyogrio
import geopandas as gpd
from shapely.geometry import box

from .abc_geo import GISDataSource
from .geo_cache import LMDBCache

class Vector(GISDataSource, LMDBCache):
    def __init__(self, filepath: str):
        self.filepath = os.path.abspath(filepath)

        LMDBCache.__init__(self, cache_name="vector_metadata")
        self._metadata = self._load_or_compute_metadata()

    @property
    def bbox(self) -> tuple:
        return tuple(self._metadata["bbox"])

    @property
    def epsg_4326_bbox(self) -> tuple:
        return tuple(self._metadata["epsg_4326_bbox"])
    
    @property
    def projection(self) -> str:
        return self._metadata["projection"]

    def _load_or_compute_metadata(self) -> dict[str, Any]:
        if not self.can_cache:
            return self._compute_metadata()

        timestamp = self._get_file_timestamp(self.filepath)

        # Fast path
        metadata = self._load_cached_metadata(
            self.filepath,
            timestamp,
        )

        if metadata is not None:
            return metadata

        # Prevent duplicate computation within process threads
        with self._CACHE_LOCK:
            metadata = self._load_cached_metadata(
                self.filepath,
                timestamp,
            )

            if metadata is not None:
                return metadata

            metadata = self._compute_metadata()

            self._save_cached_metadata(
                self.filepath,
                timestamp,
                metadata,
            )

        return metadata
    
    def _compute_metadata(self) -> dict[str, Any]:
        info = pyogrio.read_info(self.filepath, force_total_bounds=True)
        projection = info['crs']
        if projection is None and self.filepath.lower().endswith(('.parquet', '.geoparquet')):
            import pyarrow.parquet as pq

            parquet_file = pq.ParquetFile(self.filepath)

            # Retrieve the file metadata
            metadata = parquet_file.metadata
            geo_metadata = json.loads(metadata.metadata[b'geo'].decode('utf-8'))
            projection = ":".join(map(str,geo_metadata['columns']['geometry']['crs']['id'].values()))

        bbox = info['total_bounds']
        if projection is not None and projection != 'EPSG:4326':
            minx, miny, maxx, maxy = bbox
            gdf_bbox = gpd.GeoSeries([box(minx, miny, maxx, maxy)], crs=projection).to_crs(4326)
            epsg_4326_bbox = gdf_bbox.total_bounds
        else:
            epsg_4326_bbox = bbox

        return {
            "projection": projection,
            "bbox": bbox,
            "epsg_4326_bbox": epsg_4326_bbox,
        }

    def to_geopandas(self, bbox_epsg_4326: list[float] = None, columns: list[str] = None) -> gpd.GeoDataFrame:

        """
        Read any geometry file (shapefile, geojson, geoparquet) into a GeoDataFrame. Bbox in 4326.
        """
        if bbox_epsg_4326 is not None:
            # If bbox is provided, we need to make sure it's in the same CRS as the data. We can check the CRS of the data by reading just the metadata with pyogrio, and then reprojecting the bbox if necessary.
            data_crs = self.projection
            if data_crs is not None and data_crs != 'EPSG:4326':
                minx, miny, maxx, maxy = bbox_epsg_4326
                gdf_bbox = gpd.GeoSeries([box(minx, miny, maxx, maxy)], crs="EPSG:4326").to_crs(data_crs)
                bbox_epsg_4326 = gdf_bbox.total_bounds

        if self.filepath.lower().endswith(('.parquet', '.geoparquet')):
            try:
                return gpd.read_parquet(self.filepath, bbox=bbox_epsg_4326, columns=columns)
            except ValueError as e:
                if "Specifying 'bbox' not supported for this Parquet file" in str(e):
                    warnings.warn(f"Could not read {self.filepath} with bbox. Consider adding covering bbox to parquet file for faster reading.", stacklevel=2)
                    gdf = gpd.read_parquet(self.filepath, columns=columns)
                    minx, miny, maxx, maxy = bbox_epsg_4326
                    gdf = gdf.cx[minx:maxx, miny:maxy]
                    return gdf
                raise e
        
        return gpd.read_file(self.filepath, use_arrow=True, bbox=bbox_epsg_4326, columns=columns)
    
    @classmethod
    def save_any_geom(cls, gdf: gpd.GeoDataFrame, path: str, **kwargs) -> None:
        if path.lower().endswith(('.parquet', '.geoparquet')):
            gdf.to_parquet(path, **kwargs)
        else:
            gdf.to_file(path, **kwargs)

    @classmethod
    def _streamline_is_in_dem_bounds(cls, stream: str, dem_bounds: tuple[float, float, float, float]) -> bool:
        stream_bounds = Vector(stream).epsg_4326_bbox

        return cls.bounds_intersect(stream_bounds, dem_bounds)

    @classmethod
    def get_streamlines_in_extent(cls, bounds: tuple[float, float, float, float], streamlines: list[str]) -> list[str]:
        """Return the stream parquet files whose stored bounds intersect a DEM tile."""
        streamlines_to_clip = []
        for stream in streamlines:
            if cls._streamline_is_in_dem_bounds(stream, bounds):
                streamlines_to_clip.append(stream)

        return streamlines_to_clip