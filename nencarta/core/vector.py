from __future__ import annotations

import json
import os
import warnings
from pathlib import Path
from typing import Any

import pyogrio
import geopandas as gpd
from pyproj import CRS
from shapely.geometry import box

from .abc_geo import GISDataSource
from .geo_cache import LMDBCache

class Vector(GISDataSource, LMDBCache):
    def __init__(self, filepath: str, use_threads: bool = True):
        self.filepath = Path(filepath).resolve()
        self.use_threads = use_threads

        LMDBCache.__init__(self)
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

        projection = info["crs"]

        # GeoParquet may have CRS metadata even if pyogrio doesn't report it.
        if projection is None:
            if Path(self.filepath).suffix.lower() == ".parquet":
                import pyarrow.parquet as pq

                metadata = pq.ParquetFile(self.filepath).metadata.metadata or {}
                geo = metadata.get(b"geo")
                if geo:
                    crs = json.loads(geo.decode("utf-8"))["columns"]["geometry"].get("crs")
                    if crs and "id" in crs:
                        projection = ":".join(map(str, crs["id"].values()))
                    if crs and not projection:
                        projection = CRS.from_user_input(crs)
            else:
                warnings.warn(f"Could not determine CRS for {self.filepath}. Assuming EPSG:4326.", stacklevel=2)
                projection = CRS.from_epsg(4326)

        bbox = info["total_bounds"]

        epsg_4326_bbox = bbox
        if projection is not None:
            crs = CRS.from_user_input(projection)
            if crs != CRS.from_epsg(4326):
                epsg_4326_bbox = (
                    gpd.GeoSeries([box(*bbox)], crs=crs)
                    .to_crs(4326)
                    .total_bounds
                )

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

        if self.filepath.suffix.lower() in ('.parquet', '.geoparquet'):
            try:
                return gpd.read_parquet(self.filepath, bbox=bbox_epsg_4326, columns=columns, use_threads=self.use_threads)
            except ValueError as e:
                if "Specifying 'bbox' not supported for this Parquet file" in str(e):
                    # warnings.warn(f"Could not read {self.filepath} with bbox. Consider adding covering bbox to parquet file for faster reading.", stacklevel=2)
                    gdf = gpd.read_parquet(self.filepath, columns=columns, use_threads=self.use_threads)
                    minx, miny, maxx, maxy = bbox_epsg_4326
                    gdf = gdf.cx[minx:maxx, miny:maxy]
                    return gdf
                raise e
            
        if self.filepath.suffix.lower() == ".gdb":
            # Specify the layer you want to access
            layer = "geoglowsv2"
        else:
            layer = None
        
        return gpd.read_file(self.filepath, use_arrow=True, bbox=bbox_epsg_4326, columns=columns, layer=layer)
    
    @classmethod
    def save_any_geom(cls, gdf: gpd.GeoDataFrame, path: os.PathLike, **kwargs) -> None:
        cls.delete_any_geom(path)
        if Path(path).suffix.lower() in ('.parquet', '.geoparquet'):
            gdf.to_parquet(path, **kwargs)
        else:
            gdf.to_file(path, **kwargs)

    @classmethod
    def delete_any_geom(cls, path: os.PathLike) -> None:
        path = Path(path).resolve()
        if not path.exists():
            return

        if path.suffix.lower() == ".shp":
            shapefile = Path(path)
            for sidecar in shapefile.parent.glob(f"{shapefile.stem}.*"):
                sidecar.unlink()
            return

        path.unlink()

    @classmethod
    def _streamline_is_in_dem_bounds(cls, stream: str, dem_bounds: tuple[float, float, float, float]) -> bool:
        stream_bounds = Vector(stream).epsg_4326_bbox

        return cls.bounds_intersect(stream_bounds, dem_bounds)

    @classmethod
    def get_streamlines_in_extent(cls, bounds: tuple[float, float, float, float], streamlines: list[str]) -> list[str]:
        """Return the stream parquet files whose stored bounds intersect a DEM tile. Bounds should be in EPSG:4326."""
        streamlines_to_clip = []
        for stream in streamlines:
            if cls._streamline_is_in_dem_bounds(stream, bounds):
                streamlines_to_clip.append(stream)

        return streamlines_to_clip
