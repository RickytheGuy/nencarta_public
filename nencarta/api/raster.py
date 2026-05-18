from __future__ import annotations

import os
from typing import Any

from propcache import cached_property

import geopandas as gpd
from shapely.geometry import box
from osgeo import gdal, osr

from .abc_geo import GISDataSource
from .geo_cache import LMDBCache

gdal.UseExceptions()


class Raster(GISDataSource, LMDBCache):
    def __init__(self, filepath: str):
        self.filepath = os.path.abspath(filepath)

        LMDBCache.__init__(self, cache_name="raster_metadata")
        self._metadata = self._load_or_compute_metadata()

    @cached_property
    def ds(self) -> gdal.Dataset:
        return gdal.Open(self.filepath)

    @property
    def geotransform(self) -> tuple:
        return tuple(self._metadata["geotransform"])

    @property
    def projection(self) -> str:
        return self._metadata["projection"]

    @property
    def shape(self) -> tuple:
        return tuple(self._metadata["shape"])

    @property
    def resolution(self) -> tuple:
        return tuple(self._metadata["resolution"])

    @property
    def bbox(self) -> tuple:
        return tuple(self._metadata["bbox"])

    @property
    def epsg_4326_bbox(self) -> tuple:
        return tuple(self._metadata["epsg_4326_bbox"])

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
        ds = self.ds

        geotransform = ds.GetGeoTransform()
        projection = ds.GetProjection()

        shape = (
            ds.RasterXSize,
            ds.RasterYSize,
        )

        resolution = (
            abs(geotransform[1]),
            abs(geotransform[5]),
        )

        width, height = shape

        minx = geotransform[0]
        maxx = geotransform[0] + width * geotransform[1]
        miny = geotransform[3] + height * geotransform[5]
        maxy = geotransform[3]

        if miny > maxy:
            miny, maxy = maxy, miny

        bbox = (
            minx,
            miny,
            maxx,
            maxy,
        )

        srs = osr.SpatialReference(projection)
        srs.AutoIdentifyEPSG()

        if srs.IsGeographic() and srs.GetAuthorityCode(None) == "4326":
            epsg_4326_bbox = bbox
        else:
            epsg_4326_bbox = tuple(
                gpd.GeoSeries(
                    [box(*bbox)],
                    crs=projection,
                )
                .to_crs("EPSG:4326")
                .total_bounds
            )

        return {
            "geotransform": geotransform,
            "projection": projection,
            "shape": shape,
            "resolution": resolution,
            "bbox": bbox,
            "epsg_4326_bbox": epsg_4326_bbox,
        }

    @classmethod
    def get_rasters_in_extent(cls, epsg_4326_bounds: tuple, other_rasters: list[str]):
        output = set()
        for raster in other_rasters:
            raster_bounds = Raster(raster).epsg_4326_bbox
            if cls.bounds_intersect(epsg_4326_bounds, raster_bounds):
                output.add(raster)

        return list(output)