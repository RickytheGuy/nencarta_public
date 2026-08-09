# built-in imports
import io
import json
import logging
import os
import sys
import traceback
import urllib.request
from pathlib import Path

# third-party imports
from PyQt5.QtCore import QSettings, QThread, QTimer, Qt, pyqtSignal
from PyQt5.QtGui import QColor, QFont, QIcon, QImage, QPalette, QPainterPath, QPen, QPixmap, QTextCursor
from PyQt5.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QGridLayout,
    QGraphicsPathItem,
    QGraphicsPixmapItem,
    QGraphicsScene,
    QGraphicsView,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSpinBox,
    QSplitter,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

# local imports
from nencarta.core.configs import NencartaConfig
from nencarta.core.defaults import (
    DEFAULT_BATHY_ARGS,
    DEFAULT_CONFIG,
    DEFAULT_FLOODMAP_ARGS,
    FLOW_FIELD_OPTIONS,
    VALID_RETURN_PERIODS,
)
from nencarta.core.enumerations import FloodMapMode, Mapper, StreamflowSource
from nencarta.logger import LOG


settings = QSettings("NenCarta", "FloodSimulationGUI")

STREAMFLOW_SOURCE_LABELS = {
    str(StreamflowSource.GEOGLOWS): str(StreamflowSource.GEOGLOWS),
    "NWM Short Range": str(StreamflowSource.NWM_SHORT_RANGE),
    "NWM Medium Range": str(StreamflowSource.NWM_MEDIUM_RANGE),
    "NWM Long Range": str(StreamflowSource.NWM_LONG_RANGE),
}

def _streamflow_label(value):
    for label, source in STREAMFLOW_SOURCE_LABELS.items():
        if value in {label, source}:
            return label
    return str(StreamflowSource.GEOGLOWS)


class QtLogStream(io.TextIOBase):
    def __init__(self, qt_signal):
        super().__init__()
        self.qt_signal = qt_signal
        self._buffer = ""
        self._replace_next = False

    def write(self, text):
        if not text:
            return 0

        for char in text:
            if char == "\r":
                if self._buffer.strip():
                    self.qt_signal.emit(("\r" if self._replace_next else "") + self._buffer)
                self._buffer = ""
                self._replace_next = True
            elif char == "\n":
                if self._buffer.strip():
                    self.qt_signal.emit(("\r" if self._replace_next else "") + self._buffer)
                self._buffer = ""
                self._replace_next = False
            else:
                self._buffer += char

        if self._replace_next and self._buffer.strip():
            self.qt_signal.emit("\r" + self._buffer)
            self._buffer = ""
        return len(text)

    def flush(self):
        if self._buffer.strip():
            self.qt_signal.emit(("\r" if self._replace_next else "") + self._buffer)
        self._buffer = ""
        self._replace_next = False


class QtLogHandler(logging.Handler):
    def __init__(self, signal):
        super().__init__()
        self.signal = signal

    def emit(self, record):
        try:
            self.signal.emit(self.format(record))
        except Exception:
            self.handleError(record)


class MapGraphicsView(QGraphicsView):
    view_changed = pyqtSignal()

    def __init__(self):
        super().__init__()
        self.setRenderHints(self.renderHints())
        self.setDragMode(QGraphicsView.ScrollHandDrag)
        self.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.AnchorViewCenter)

    def wheelEvent(self, event):
        factor = 1.25 if event.angleDelta().y() > 0 else 0.8
        self.scale(factor, factor)
        self.view_changed.emit()

    def mouseReleaseEvent(self, event):
        super().mouseReleaseEvent(event)
        self.view_changed.emit()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.view_changed.emit()


class MapLayerLoader:
    RASTER_SUFFIXES = {".tif", ".tiff", ".vrt", ".img"}
    VECTOR_SUFFIXES = {".gpkg", ".shp", ".gdb", ".parquet", ".geoparquet", ".geojson", ".json"}
    COLORS = ("#2563eb", "#16a34a", "#dc2626", "#9333ea", "#ea580c", "#0891b2", "#4f46e5")
    EXCLUDED_OUTPUT_PATTERNS = (
        "_water_mask",
        "_strm_raster",
        "_strm_raster_clean",
        "_matched.tif",
        "_wtbx_derived.tif",
    )

    @classmethod
    def layer_kind(cls, path):
        suffix = Path(path).suffix.lower()
        if suffix in cls.RASTER_SUFFIXES:
            return "raster"
        if suffix in cls.VECTOR_SUFFIXES:
            return "vector"
        return "file"

    @classmethod
    def discover_rasters(cls, directory):
        directory = Path(directory)
        if not directory.is_dir():
            return []
        paths = []
        for suffix in cls.RASTER_SUFFIXES:
            paths.extend(directory.glob(f"*{suffix}"))
        return sorted(paths)

    @classmethod
    def discover_outputs(cls, output_dir, watershed_name):
        root = Path(output_dir) / watershed_name
        if not root.is_dir():
            return []
        suffixes = cls.RASTER_SUFFIXES | cls.VECTOR_SUFFIXES
        return sorted(
            path
            for path in root.rglob("*")
            if path.suffix.lower() in suffixes and not cls.skip_output_layer(path)
        )

    @classmethod
    def skip_output_layer(cls, path):
        name = Path(path).name.lower()
        return any(pattern in name for pattern in cls.EXCLUDED_OUTPUT_PATTERNS)

    @classmethod
    def google_hybrid_tile_specs(cls, bbox, max_tiles=64):
        min_lon, min_lat, max_lon, max_lat = bbox
        if min_lon >= max_lon or min_lat >= max_lat:
            return []

        zoom = cls.choose_tile_zoom(max_lon - min_lon, max_lat - min_lat)
        x_min, y_max = cls.lonlat_to_tile(min_lon, min_lat, zoom)
        x_max, y_min = cls.lonlat_to_tile(max_lon, max_lat, zoom)
        x_min, x_max = sorted((x_min, x_max))
        y_min, y_max = sorted((y_min, y_max))
        center_x, center_y = cls.lonlat_to_tile((min_lon + max_lon) / 2, (min_lat + max_lat) / 2, zoom)

        tiles = []
        for x in range(x_min, x_max + 1):
            for y in range(y_min, y_max + 1):
                tiles.append((zoom, x, y))
        tiles.sort(key=lambda tile: (tile[1] - center_x) ** 2 + (tile[2] - center_y) ** 2)
        return tiles[:max_tiles]

    @staticmethod
    def choose_tile_zoom(width_degrees, height_degrees):
        span = max(width_degrees, height_degrees)
        if span > 20:
            return 5
        if span > 8:
            return 6
        if span > 3:
            return 8
        if span > 1:
            return 10
        if span > 0.25:
            return 12
        if span > 0.05:
            return 14
        return 16

    @staticmethod
    def lonlat_to_tile(lon, lat, zoom):
        import math

        lat = max(min(lat, 85.05112878), -85.05112878)
        n = 2 ** zoom
        x = int((lon + 180.0) / 360.0 * n)
        lat_rad = math.radians(lat)
        y = int((1.0 - math.asinh(math.tan(lat_rad)) / math.pi) / 2.0 * n)
        return max(0, min(n - 1, x)), max(0, min(n - 1, y))

    @staticmethod
    def tile_bounds(x, y, zoom):
        import math

        n = 2 ** zoom
        lon_min = x / n * 360.0 - 180.0
        lon_max = (x + 1) / n * 360.0 - 180.0
        lat_max = math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * y / n))))
        lat_min = math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * (y + 1) / n))))
        return lon_min, lat_min, lon_max, lat_max

    @staticmethod
    def google_tile_url(x, y, zoom):
        return f"https://mt1.google.com/vt/lyrs=y&x={x}&y={y}&z={zoom}"

    @classmethod
    def make_google_tile_item_from_data(cls, x, y, zoom, data):
        pixmap = QPixmap()
        if not pixmap.loadFromData(data):
            return None

        item = QGraphicsPixmapItem(pixmap)
        minx, miny, maxx, maxy = cls.tile_bounds(x, y, zoom)
        item.setPos(minx, -maxy)
        item.setTransform(item.transform().scale((maxx - minx) / pixmap.width(), (maxy - miny) / pixmap.height()), True)
        item.setZValue(-100)
        return item

    @classmethod
    def make_google_tile_item(cls, x, y, zoom):
        url = cls.google_tile_url(x, y, zoom)
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "NenCarta GUI"})
            with urllib.request.urlopen(request, timeout=2) as response:
                data = response.read()
        except Exception:
            return None
        return cls.make_google_tile_item_from_data(x, y, zoom, data)

    @classmethod
    def make_layer(cls, name, path, group, index):
        path = Path(path)
        kind = cls.layer_kind(path)
        color = cls.COLORS[index % len(cls.COLORS)]
        layer = {
            "name": name,
            "path": path,
            "group": group,
            "kind": kind,
            "color": color,
            "item": None,
            "stats": "",
            "error": None,
        }
        try:
            if not path.exists():
                layer["error"] = "File does not exist."
                layer["stats"] = cls.file_stats(path, missing=True)
            elif kind == "raster":
                layer["item"], layer["stats"] = cls.make_raster_item(path)
            elif kind == "vector":
                layer["item"], layer["stats"] = cls.make_vector_item(path, color)
            else:
                layer["stats"] = cls.file_stats(path)
        except Exception as exc:
            layer["error"] = str(exc)
            layer["stats"] = f"{path}\nError: {exc}"
        return layer

    @classmethod
    def make_raster_item(cls, path):
        import numpy as np
        from osgeo import gdal
        from pyproj import CRS

        ds = gdal.Open(str(path))
        if ds is None:
            raise ValueError("Could not open raster.")

        band = ds.GetRasterBand(1)
        x_size = ds.RasterXSize
        y_size = ds.RasterYSize
        scale = max(x_size / 900, y_size / 900, 1)
        out_x = max(1, int(x_size / scale))
        out_y = max(1, int(y_size / scale))
        array = band.ReadAsArray(buf_xsize=out_x, buf_ysize=out_y).astype("float64")
        nodata = band.GetNoDataValue()
        projection = ds.GetProjection()
        geotransform = ds.GetGeoTransform()

        bbox = cls.raster_bbox(geotransform, x_size, y_size)
        epsg_bbox = cls.to_epsg4326_bbox(bbox, projection)

        valid = np.isfinite(array)
        if nodata is not None:
            valid &= array != nodata

        color_table = band.GetRasterColorTable()
        color_table_count = color_table.GetCount() if color_table is not None else 0
        if color_table_count:
            rgba = np.zeros((array.shape[0], array.shape[1], 4), dtype="uint8")
            indices = np.rint(array).astype("int64")
            for value in np.unique(indices[valid]):
                if 0 <= value < color_table_count:
                    red, green, blue, alpha = color_table.GetColorEntry(int(value))
                    rgba[indices == value] = [red, green, blue, alpha]
            rgba[..., 3] = np.where(valid, rgba[..., 3], 0).astype("uint8")
        else:
            if valid.any():
                low, high = np.nanpercentile(array[valid], [2, 98])
                if low == high:
                    high = low + 1
                normalized = np.clip((array - low) / (high - low), 0, 1)
                gray = (normalized * 255).astype("uint8")
            else:
                gray = np.zeros(array.shape, dtype="uint8")

            rgba = np.zeros((array.shape[0], array.shape[1], 4), dtype="uint8")
            rgba[..., 0] = gray
            rgba[..., 1] = gray
            rgba[..., 2] = gray
            rgba[..., 3] = np.where(valid, 210, 0).astype("uint8")

        image = QImage(rgba.data, rgba.shape[1], rgba.shape[0], QImage.Format_RGBA8888).copy()
        pixmap = QPixmap.fromImage(image)
        item = QGraphicsPixmapItem(pixmap)
        minx, miny, maxx, maxy = epsg_bbox
        item.setPos(minx, -maxy)
        item.setScale(1)
        item.setTransformOriginPoint(0, 0)
        item.setTransform(item.transform().scale((maxx - minx) / out_x, (maxy - miny) / out_y), True)
        item.setZValue(0)

        crs_text = CRS.from_user_input(projection).to_string() if projection else "Unknown"
        stats = "\n".join(
            [
                f"{path}",
                "Type: Raster",
                f"Projection: {crs_text}",
                f"Size: {x_size} x {y_size}",
                f"Resolution: {abs(geotransform[1]):g}, {abs(geotransform[5]):g}",
                f"BBox EPSG:4326: {cls.format_bbox(epsg_bbox)}",
                f"NoData: {nodata}",
                f"Color table entries: {color_table_count}",
            ]
        )
        ds = None
        return item, stats

    @classmethod
    def make_vector_item(cls, path, color):
        import geopandas as gpd
        import pyogrio
        from pyproj import CRS

        info = pyogrio.read_info(path, force_total_bounds=True)
        if Path(path).suffix.lower() in {".parquet", ".geoparquet"}:
            gdf = gpd.read_parquet(path)
            if len(gdf) > 5000:
                gdf = gdf.head(5000)
        else:
            kwargs = {"max_features": 5000}
            if Path(path).suffix.lower() == ".gdb":
                kwargs["layer"] = "geoglowsv2"
            gdf = pyogrio.read_dataframe(path, **kwargs)

        if gdf.empty:
            raise ValueError("Vector file has no features to draw.")
        if gdf.crs is not None:
            gdf = gdf.to_crs("EPSG:4326")

        painter_path = QPainterPath()
        for geom in gdf.geometry:
            cls.add_geometry_to_path(painter_path, geom)

        item = QGraphicsPathItem(painter_path)
        item.setPen(QPen(QColor(color), 0))
        item.setZValue(10)

        bounds = tuple(gdf.total_bounds)
        crs_text = CRS.from_user_input(info["crs"]).to_string() if info.get("crs") is not None else "Unknown"
        feature_count = info.get("features")
        stats = "\n".join(
            [
                f"{path}",
                "Type: Vector",
                f"Projection: {crs_text}",
                f"Geometry: {info.get('geometry_type')}",
                f"Features: {feature_count}",
                f"Displayed: {len(gdf)}",
                f"BBox EPSG:4326: {cls.format_bbox(bounds)}",
            ]
        )
        return item, stats

    @staticmethod
    def add_geometry_to_path(path, geom):
        if geom is None or geom.is_empty:
            return
        geom_type = geom.geom_type
        if geom_type == "LineString":
            MapLayerLoader.add_line(path, geom.coords)
        elif geom_type == "MultiLineString":
            for line in geom.geoms:
                MapLayerLoader.add_line(path, line.coords)
        elif geom_type == "Polygon":
            MapLayerLoader.add_polygon(path, geom)
        elif geom_type == "MultiPolygon":
            for polygon in geom.geoms:
                MapLayerLoader.add_polygon(path, polygon)
        elif geom_type == "Point":
            path.addEllipse(geom.x - 0.0005, -geom.y - 0.0005, 0.001, 0.001)
        elif geom_type == "MultiPoint":
            for point in geom.geoms:
                path.addEllipse(point.x - 0.0005, -point.y - 0.0005, 0.001, 0.001)
        elif hasattr(geom, "geoms"):
            for part in geom.geoms:
                MapLayerLoader.add_geometry_to_path(path, part)

    @staticmethod
    def add_line(path, coords):
        coords = list(coords)
        if not coords:
            return
        path.moveTo(coords[0][0], -coords[0][1])
        for x, y, *_ in coords[1:]:
            path.lineTo(x, -y)

    @staticmethod
    def add_polygon(path, polygon):
        MapLayerLoader.add_line(path, polygon.exterior.coords)
        for interior in polygon.interiors:
            MapLayerLoader.add_line(path, interior.coords)

    @staticmethod
    def raster_bbox(geotransform, width, height):
        xs = [
            geotransform[0],
            geotransform[0] + width * geotransform[1],
            geotransform[0] + height * geotransform[2],
            geotransform[0] + width * geotransform[1] + height * geotransform[2],
        ]
        ys = [
            geotransform[3],
            geotransform[3] + width * geotransform[4],
            geotransform[3] + height * geotransform[5],
            geotransform[3] + width * geotransform[4] + height * geotransform[5],
        ]
        return min(xs), min(ys), max(xs), max(ys)

    @staticmethod
    def to_epsg4326_bbox(bbox, projection):
        if not projection:
            return bbox
        import geopandas as gpd
        from pyproj import CRS
        from shapely.geometry import box

        crs = CRS.from_user_input(projection)
        if crs == CRS.from_epsg(4326):
            return bbox
        return tuple(gpd.GeoSeries([box(*bbox)], crs=crs).to_crs("EPSG:4326").total_bounds)

    @staticmethod
    def format_bbox(bbox):
        return ", ".join(f"{value:.6g}" for value in bbox)

    @staticmethod
    def file_stats(path, missing=False):
        size = "missing" if missing else f"{Path(path).stat().st_size / 1024:.1f} KB"
        return "\n".join([f"{path}", "Type: File", f"Size: {size}"])


class BasemapTileWorker(QThread):
    tile_signal = pyqtSignal(int, int, int, int, object)
    finished_signal = pyqtSignal(int, int)

    def __init__(self, request_id, tile_specs, max_workers=8):
        super().__init__()
        self.request_id = request_id
        self.tile_specs = list(tile_specs)
        self.max_workers = max(1, int(max_workers))

    @staticmethod
    def fetch_tile(tile_spec):
        zoom, x, y = tile_spec
        url = MapLayerLoader.google_tile_url(x, y, zoom)
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "NenCarta GUI"})
            with urllib.request.urlopen(request, timeout=1) as response:
                data = response.read()
        except Exception:
            return None
        return zoom, x, y, data

    def run(self):
        import concurrent.futures

        loaded = 0
        tile_iter = iter(self.tile_specs)
        max_workers = min(self.max_workers, len(self.tile_specs))
        executor = concurrent.futures.ThreadPoolExecutor(max_workers=max_workers)
        pending = set()

        def submit_next():
            if self.isInterruptionRequested():
                return False
            try:
                tile_spec = next(tile_iter)
            except StopIteration:
                return False
            pending.add(executor.submit(self.fetch_tile, tile_spec))
            return True

        try:
            for _ in range(max_workers):
                submit_next()

            while pending and not self.isInterruptionRequested():
                done, pending = concurrent.futures.wait(
                    pending,
                    timeout=0.1,
                    return_when=concurrent.futures.FIRST_COMPLETED,
                )
                for future in done:
                    result = future.result()
                    if result is not None and not self.isInterruptionRequested():
                        zoom, x, y, data = result
                        loaded += 1
                        self.tile_signal.emit(self.request_id, zoom, x, y, data)
                    submit_next()
        finally:
            for future in pending:
                future.cancel()
            executor.shutdown(wait=False, cancel_futures=True)
        self.finished_signal.emit(self.request_id, loaded)


class WorkerThread(QThread):
    finished_signal = pyqtSignal(str)
    log_signal = pyqtSignal(str)
    error_signal = pyqtSignal(str)

    def __init__(self, params: dict):
        super().__init__()
        self.params = dict(params)

    @staticmethod
    def _format_exception(e: Exception, context: str = ""):
        exc_type = type(e).__name__
        exc_msg = str(e)
        location = ""
        code_line = ""
        try:
            frames = traceback.extract_tb(e.__traceback__)
            if frames:
                last = frames[-1]
                location = f"{os.path.basename(last.filename)}:{last.lineno} in {last.name}()"
                code_line = (last.line or "").strip()
        except Exception:
            pass

        summary_parts = [context] if context else []
        summary_parts.append(f"{exc_type}: {exc_msg}".strip())
        if location:
            summary_parts.append(f"Location: {location}")
        if code_line:
            summary_parts.append(f"Code: {code_line}")
        return "\n".join(summary_parts), traceback.format_exc()

    def run(self):
        qt_handler = QtLogHandler(self.log_signal)
        qt_handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        LOG.addHandler(qt_handler)
        LOG.setLevel(logging.INFO)

        old_stdout, old_stderr = sys.stdout, sys.stderr
        stream = QtLogStream(self.log_signal)

        try:
            sys.stdout = stream
            sys.stderr = stream
            from nencarta.main import process_watershed

            LOG.info("Running nencarta.process_watershed()")
            LOG.info(json.dumps(self.params, indent=2, default=str))
            process_watershed(dict(self.params))
            self.finished_signal.emit(self.params["name"])
        except Exception as e:
            summary, details = self._format_exception(e, context="Error running simulation")
            LOG.error(summary.replace("\n", " | "))
            self.error_signal.emit(summary + "\n\nDETAILS:\n" + details)
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr
            LOG.removeHandler(qt_handler)


class PathPicker(QWidget):
    def __init__(self, mode="file", default_path="", dialog_title="Select Path", file_filter="All Files (*)"):
        super().__init__()
        self.mode = mode
        self.dialog_title = dialog_title
        self.file_filter = file_filter
        self.line_edit = QLineEdit(str(default_path or ""))
        self.browse_btn = QPushButton("Browse...")
        self.browse_btn.clicked.connect(self._browse)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        layout.addWidget(self.line_edit)
        layout.addWidget(self.browse_btn)

    def _browse(self):
        current = self.line_edit.text().strip()
        start_dir = _dialog_start_dir(current)
        if self.mode == "directory":
            value = QFileDialog.getExistingDirectory(
                self,
                self.dialog_title,
                start_dir,
                _dialog_options(QFileDialog.ShowDirsOnly, QFileDialog.DontResolveSymlinks),
            )
        else:
            value, _ = QFileDialog.getOpenFileName(
                self,
                self.dialog_title,
                start_dir,
                self.file_filter,
                options=_dialog_options(QFileDialog.DontResolveSymlinks),
            )
        if value:
            _remember_dialog_dir(value)
            self.line_edit.setText(value)

    def text(self):
        return self.line_edit.text()


class MultiFilePicker(QWidget):
    changed_signal = pyqtSignal()

    def __init__(self, default_paths=None, dialog_title="Select Files", file_filter="All Files (*)"):
        super().__init__()
        self.dialog_title = dialog_title
        self.file_filter = file_filter
        self.list_widget = QListWidget()
        self.list_widget.setSelectionMode(QListWidget.ExtendedSelection)
        self.list_widget.setMinimumHeight(96)

        add_btn = QPushButton("Add files...")
        remove_btn = QPushButton("Remove")
        clear_btn = QPushButton("Clear")
        add_btn.clicked.connect(self._browse)
        remove_btn.clicked.connect(self.remove_selected)
        clear_btn.clicked.connect(self.clear)

        buttons = QHBoxLayout()
        buttons.addWidget(add_btn)
        buttons.addWidget(remove_btn)
        buttons.addWidget(clear_btn)
        buttons.addStretch(1)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        layout.addWidget(self.list_widget)
        layout.addLayout(buttons)

        self.set_values(default_paths or [])

    def _browse(self):
        current = self.values()[-1] if self.values() else ""
        files, _ = QFileDialog.getOpenFileNames(
            self,
            self.dialog_title,
            _dialog_start_dir(current),
            self.file_filter,
            options=_dialog_options(QFileDialog.DontResolveSymlinks),
        )
        if files:
            _remember_dialog_dir(files[0])
        for file_path in files:
            self.add_value(file_path)

    def add_value(self, value):
        value = str(value).strip()
        if value and value not in self.values():
            self.list_widget.addItem(value)
            self.changed_signal.emit()

    def remove_selected(self):
        for item in self.list_widget.selectedItems():
            self.list_widget.takeItem(self.list_widget.row(item))
        self.changed_signal.emit()

    def clear(self):
        self.list_widget.clear()
        self.changed_signal.emit()

    def set_values(self, values):
        old_values = self.values()
        self.clear()
        if isinstance(values, str):
            values = _split_list(values)
        for value in values:
            self.add_value(value)
        if self.values() != old_values:
            self.changed_signal.emit()

    def values(self):
        return [self.list_widget.item(row).text() for row in range(self.list_widget.count())]


class ReturnPeriodPicker(QWidget):
    def __init__(self, selected=None):
        super().__init__()
        selected = {int(value) for value in (selected or [])}
        self.checkboxes = {}
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(12)
        for period in VALID_RETURN_PERIODS:
            checkbox = QCheckBox(str(period))
            checkbox.setChecked(period in selected)
            self.checkboxes[period] = checkbox
            layout.addWidget(checkbox)
        layout.addStretch(1)

    def values(self):
        return [period for period, checkbox in self.checkboxes.items() if checkbox.isChecked()]

    def set_values(self, values):
        values = {int(value) for value in (values or [])}
        for period, checkbox in self.checkboxes.items():
            checkbox.setChecked(period in values)


class DictTable(QWidget):
    def __init__(self, data=None, defaults=None):
        super().__init__()
        self.defaults = dict(defaults or {})
        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["Parameter", "Value"])
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.table.verticalHeader().setVisible(False)
        self.table.setAlternatingRowColors(True)
        self.table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        add_btn = QPushButton("Add argument")
        remove_btn = QPushButton("Remove")
        reset_btn = QPushButton("Reset defaults")
        add_btn.clicked.connect(self.add_row)
        remove_btn.clicked.connect(self.remove_selected_rows)
        reset_btn.clicked.connect(self.reset_defaults)

        buttons = QHBoxLayout()
        buttons.addWidget(add_btn)
        buttons.addWidget(remove_btn)
        buttons.addWidget(reset_btn)
        buttons.addStretch(1)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.table)
        layout.addLayout(buttons)

        self.set_dict(data or {})

    def add_row(self, key="", value=""):
        row = self.table.rowCount()
        self.table.insertRow(row)
        self.table.setItem(row, 0, QTableWidgetItem(str(key)))
        self.table.setItem(row, 1, QTableWidgetItem(str(value)))

    def remove_selected_rows(self):
        rows = sorted({index.row() for index in self.table.selectedIndexes()}, reverse=True)
        for row in rows:
            self.table.removeRow(row)

    def set_dict(self, data: dict):
        self.table.setRowCount(0)
        for key, value in data.items():
            self.add_row(key, _format_setting_value(value))

    def reset_defaults(self):
        self.set_dict(self.defaults)

    def to_dict(self) -> dict:
        result = {}
        for row in range(self.table.rowCount()):
            key_item = self.table.item(row, 0)
            value_item = self.table.item(row, 1)
            if key_item and key_item.text().strip():
                result[key_item.text().strip()] = _parse_scalar(value_item.text() if value_item else "")
        return result


def _parse_scalar(value):
    value = str(value).strip()
    if value == "":
        return ""
    lowered = value.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    if lowered in {"none", "null"}:
        return None
    try:
        return json.loads(value)
    except Exception:
        pass
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def _format_setting_value(value):
    if isinstance(value, (dict, list, bool, int, float)):
        return json.dumps(value)
    return "" if value is None else str(value)


def _setting_text(key, default="", multiline=False):
    value = settings.value(key, default)
    if value is None:
        return ""
    if isinstance(value, str):
        try:
            parsed = json.loads(value)
        except (TypeError, json.JSONDecodeError):
            return value
        if isinstance(parsed, list):
            separator = "\n" if multiline else ", "
            return separator.join(str(item) for item in parsed)
        return str(parsed)
    if isinstance(value, list):
        separator = "\n" if multiline else ", "
        return separator.join(str(item) for item in value)
    return str(value)


def _setting_list(key, default=None):
    value = settings.value(key, default if default is not None else [])
    if value is None:
        return []
    if isinstance(value, str):
        try:
            parsed = json.loads(value)
        except (TypeError, json.JSONDecodeError):
            return _split_list(value)
        if isinstance(parsed, list):
            return parsed
        return [str(parsed)]
    if isinstance(value, list):
        return value
    return [value]


def _default(key, fallback=""):
    return DEFAULT_CONFIG.get(key, fallback)


def _split_list(value):
    text = str(value or "").replace(",", "\n")
    return [part.strip() for part in text.splitlines() if part.strip()]


def _parse_float_list(value):
    return [float(part) for part in _split_list(value)]


def _parse_int_list(value):
    return [int(float(part)) for part in _split_list(value)]


def _parse_bbox(value):
    parts = _parse_float_list(value)
    if not parts:
        return None
    if len(parts) != 4:
        raise ValueError("bbox must contain four values: minx, miny, maxx, maxy.")
    return parts


def _is_network_path(path_text):
    text = str(path_text or "").strip()
    return text.startswith("\\\\") or text.startswith("//")


def _local_dialog_dir(path_text):
    text = str(path_text or "").strip()
    if not text or _is_network_path(text):
        return None
    try:
        path = Path(text).expanduser()
        if path.is_dir():
            return str(path)
        if path.is_file():
            return str(path.parent)
        if path.parent.is_dir():
            return str(path.parent)
    except (OSError, RuntimeError, ValueError):
        return None
    return None


def _dialog_start_dir(path_text=""):
    for candidate in (
        path_text,
        settings.value("last_file_dialog_dir", ""),
        str(Path.home()),
    ):
        directory = _local_dialog_dir(candidate)
        if directory:
            return directory
    return ""


def _remember_dialog_dir(path_text):
    if not path_text or _is_network_path(path_text):
        return
    try:
        path = Path(path_text).expanduser()
        directory = path if path.is_dir() else path.parent
        if directory.is_dir():
            settings.setValue("last_file_dialog_dir", str(directory))
    except (OSError, RuntimeError, ValueError):
        return


def _dialog_options(*options):
    combined = QFileDialog.DontUseNativeDialog
    for option in options:
        combined |= option
    return combined


class FloodSimulationGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("NenCarta")
        icon_path = Path(__file__).with_name("images") / "header_logo.png"
        if icon_path.exists():
            self.setWindowIcon(QIcon(str(icon_path)))

        self.resize(1360, 860)
        self.input_fields = {}
        self.field_parsers = {}
        self.worker_thread = None
        self._last_log_message_was_progress = False
        self._init_ui()

    def _init_ui(self):
        self.central_widget = QWidget()
        self.central_widget.setObjectName("centralWidget")
        self.setCentralWidget(self.central_widget)

        root_layout = QVBoxLayout(self.central_widget)
        root_layout.setContentsMargins(18, 16, 18, 16)
        root_layout.setSpacing(14)

        header = QWidget()
        header_layout = QHBoxLayout(header)
        header_layout.setContentsMargins(0, 0, 0, 0)

        title_block = QVBoxLayout()
        title = QLabel("NenCarta")
        title.setObjectName("AppTitle")
        subtitle = QLabel("Flood simulation setup")
        subtitle.setObjectName("AppSubtitle")
        title_block.addWidget(title)
        title_block.addWidget(subtitle)

        self.run_button = QPushButton("Start Simulation")
        self.run_button.setObjectName("PrimaryButton")
        self.run_button.clicked.connect(self.start_simulation)
        self.run_button.setMinimumHeight(40)

        header_layout.addLayout(title_block)
        header_layout.addStretch(1)
        header_layout.addWidget(self.run_button)
        root_layout.addWidget(header)

        splitter = QSplitter(Qt.Horizontal)
        root_layout.addWidget(splitter, 1)

        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)
        self.tabs.setMinimumWidth(560)
        splitter.addWidget(self.tabs)

        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(12, 0, 0, 0)
        right_layout.setSpacing(12)
        splitter.addWidget(right_panel)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)

        self.side_tabs = QTabWidget()
        right_layout.addWidget(self.side_tabs, 1)

        log_tab = QWidget()
        log_layout = QVBoxLayout(log_tab)
        log_layout.setContentsMargins(0, 0, 0, 0)
        log_group = QGroupBox("Simulation Log")
        log_group_layout = QVBoxLayout(log_group)
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setFont(QFont("Consolas", 9))
        log_group_layout.addWidget(self.log_text)
        log_layout.addWidget(log_group)
        self.side_tabs.addTab(log_tab, "Log")

        self._init_map_panel()

        self._add_input_tabs()
        self._wire_visibility()
        self._wire_map_refresh()
        self.log_text.setText("[INFO] Configure inputs, then start the simulation.\n")

    def _init_map_panel(self):
        map_tab = QWidget()
        map_layout = QVBoxLayout(map_tab)
        map_layout.setContentsMargins(0, 0, 0, 0)
        map_layout.setSpacing(10)

        toolbar = QHBoxLayout()
        refresh_inputs = QPushButton("Refresh Inputs")
        refresh_outputs = QPushButton("Load Outputs")
        zoom_layers = QPushButton("Zoom to Layers")
        clear_layers = QPushButton("Clear")
        self.google_hybrid_basemap = QCheckBox("Google Hybrid")
        refresh_inputs.clicked.connect(self.refresh_input_layers)
        refresh_outputs.clicked.connect(self.load_output_layers)
        zoom_layers.clicked.connect(self.zoom_to_map_layers)
        clear_layers.clicked.connect(self.clear_map_layers)
        self.google_hybrid_basemap.toggled.connect(self.render_map_layers)
        toolbar.addWidget(refresh_inputs)
        toolbar.addWidget(refresh_outputs)
        toolbar.addWidget(zoom_layers)
        toolbar.addWidget(clear_layers)
        toolbar.addWidget(self.google_hybrid_basemap)
        toolbar.addStretch(1)
        map_layout.addLayout(toolbar)

        map_splitter = QSplitter(Qt.Horizontal)
        self.map_scene = QGraphicsScene(self)
        self.map_view = MapGraphicsView()
        self.map_view.setScene(self.map_scene)
        self.map_view.setMinimumHeight(360)
        map_splitter.addWidget(self.map_view)

        layer_panel = QWidget()
        layer_panel.setMinimumWidth(300)
        layer_layout = QVBoxLayout(layer_panel)
        layer_layout.setContentsMargins(0, 0, 0, 0)
        self.layer_controls = QWidget()
        self.layer_controls_layout = QVBoxLayout(self.layer_controls)
        self.layer_controls_layout.setContentsMargins(8, 8, 8, 8)
        self.layer_controls_layout.setSpacing(6)
        self.layer_controls_layout.addStretch(1)
        layer_scroll = QScrollArea()
        layer_scroll.setWidgetResizable(True)
        layer_scroll.setWidget(self.layer_controls)
        layer_layout.addWidget(QLabel("Layers"))
        layer_layout.addWidget(layer_scroll, 1)
        map_splitter.addWidget(layer_panel)
        map_splitter.setStretchFactor(0, 4)
        map_splitter.setStretchFactor(1, 1)
        map_layout.addWidget(map_splitter, 1)

        self.side_tabs.addTab(map_tab, "Map")
        self.layer_groups = {"Inputs": [], "Outputs": []}
        self.layer_visibility = {}
        self.basemap_items = []
        self.basemap_tile_cache = {}
        self.basemap_workers = []
        self.basemap_request_id = 0
        self.current_layer_rect = None
        self._basemap_warning_shown = False
        self.map_refresh_timer = QTimer(self)
        self.map_refresh_timer.setSingleShot(True)
        self.map_refresh_timer.timeout.connect(self.refresh_input_layers)
        self.basemap_refresh_timer = QTimer(self)
        self.basemap_refresh_timer.setSingleShot(True)
        self.basemap_refresh_timer.timeout.connect(self.refresh_basemap_for_view)
        self.map_view.view_changed.connect(self.schedule_basemap_refresh)
        self.map_view.horizontalScrollBar().valueChanged.connect(self.schedule_basemap_refresh)
        self.map_view.verticalScrollBar().valueChanged.connect(self.schedule_basemap_refresh)

    def _make_tab(self, title):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(12)
        layout.addStretch(1)
        scroll.setWidget(content)
        self.tabs.addTab(scroll, title)
        return layout

    def _make_group(self, parent_layout, title):
        group = QGroupBox(title)
        layout = QGridLayout(group)
        layout.setContentsMargins(14, 18, 14, 14)
        layout.setHorizontalSpacing(12)
        layout.setVerticalSpacing(10)
        layout.setColumnStretch(0, 0)
        layout.setColumnStretch(1, 1)
        parent_layout.insertWidget(parent_layout.count() - 1, group)
        return layout

    def _add_row(self, layout, label, widget, key, parser="text", tooltip=None):
        row = layout.rowCount()
        label_widget = QLabel(label)
        label_widget.setWordWrap(True)
        if tooltip:
            label_widget.setToolTip(tooltip)
            widget.setToolTip(tooltip)
        layout.addWidget(label_widget, row, 0)
        layout.addWidget(widget, row, 1)
        self.input_fields[key] = widget
        self.field_parsers[key] = parser
        return widget

    def _add_checkbox(self, layout, key, label, default=False, tooltip=None):
        widget = QCheckBox(label)
        widget.setChecked(settings.value(key, default, type=bool))
        if tooltip:
            widget.setToolTip(tooltip)
        row = layout.rowCount()
        layout.addWidget(widget, row, 0, 1, 2)
        self.input_fields[key] = widget
        self.field_parsers[key] = "bool"
        return widget

    def _add_line(self, layout, key, label, default="", parser="text", placeholder="", tooltip=None):
        value = _setting_text(key, default)
        widget = QLineEdit(value)
        widget.setPlaceholderText(placeholder)
        return self._add_row(layout, label, widget, key, parser, tooltip)

    def _add_multiline(self, layout, key, label, default="", parser="str_list", tooltip=None):
        value = _setting_text(key, default, multiline=True)
        widget = QPlainTextEdit(value)
        widget.setMinimumHeight(72)
        return self._add_row(layout, label, widget, key, parser, tooltip)

    def _add_combo(self, layout, key, label, items, default=None, editable=False, parser="text", tooltip=None):
        saved = settings.value(key, default if default is not None else "")
        widget = QComboBox()
        widget.setEditable(editable)
        widget.addItems([str(item) for item in items])
        if saved is not None:
            widget.setCurrentText(str(saved))
        return self._add_row(layout, label, widget, key, parser, tooltip)

    def _add_spin(self, layout, key, label, default=0, minimum=0, maximum=999999, tooltip=None):
        widget = QSpinBox()
        widget.setRange(minimum, maximum)
        widget.setValue(settings.value(key, default, type=int))
        return self._add_row(layout, label, widget, key, "int", tooltip)

    def _add_path(self, layout, key, label, default="", mode="file", file_filter="All Files (*)", tooltip=None):
        widget = PathPicker(mode=mode, default_path=_setting_text(key, default), dialog_title=f"Select {label}", file_filter=file_filter)
        return self._add_row(layout, label, widget, key, "path", tooltip)

    def _add_multi_files(self, layout, key, label, default=None, file_filter="All Files (*)", tooltip=None):
        widget = MultiFilePicker(
            _setting_list(key, default),
            dialog_title=f"Select {label}",
            file_filter=file_filter,
        )
        return self._add_row(layout, label, widget, key, "str_list", tooltip)

    def _add_return_periods(self, layout, key, label, default=None, tooltip=None):
        widget = ReturnPeriodPicker(_setting_list(key, default))
        return self._add_row(layout, label, widget, key, "int_list", tooltip)

    def _add_dict(self, layout, key, label, default):
        saved = settings.value(key)
        data = default
        if saved:
            try:
                data = json.loads(saved)
            except (TypeError, json.JSONDecodeError):
                data = default
        widget = DictTable(data, defaults=default)
        widget.setMinimumHeight(280)
        return self._add_row(layout, label, widget, key, "dict")

    def _add_input_tabs(self):
        inputs_tab = self._make_tab("Inputs")
        run_tab = self._make_tab("Run")
        hydro_tab = self._make_tab("Hydro")
        bathy_tab = self._make_tab("Bathymetry")
        outputs_tab = self._make_tab("Outputs")
        advanced_tab = self._make_tab("Advanced")

        source_files = "Raster/Vector Files (*.tif *.tiff *.vrt *.gpkg *.shp *.gdb *.parquet *.geoparquet);;All Files (*)"
        vector_files = "Vector Files (*.gpkg *.shp *.gdb *.parquet *.geoparquet);;All Files (*)"

        group = self._make_group(inputs_tab, "Required")
        self._add_line(group, "watershed_name", "Watershed name", settings.value("name", "example"), placeholder="example")
        self._add_path(group, "output_dir", "Output directory", "data/output", mode="directory")

        group = self._make_group(inputs_tab, "DEM Source: Directory or File")
        self._add_path(group, "dem_dir", "DEM directory", "", mode="directory")
        self._add_path(group, "dem", "DEM file", "", file_filter=source_files)

        group = self._make_group(inputs_tab, "Generated DEM")
        self._add_multi_files(group, "source_dems", "Source DEM files", _default("source_dems"), source_files)
        self._add_line(group, "bbox", "Bounding box", "", parser="bbox", placeholder="minx, miny, maxx, maxy")
        self._add_checkbox(group, "buffer", "Buffer DEM from source DEMs", _default("buffer"))
        self._add_line(group, "buffer_distance", "Buffer distance", _default("buffer_distance"), parser="float")
        self._add_checkbox(group, "use_vrt", "Create assigned DEM as VRT", _default("use_vrt"))
        self._add_checkbox(group, "use_warning_flags_to_download_dem", "Use warning flags to download DEM", _default("use_warning_flags_to_download_dem"))
        self._add_line(group, "dem_filter", "DEM filter", _default("dem_filter"), placeholder="*.tif")

        group = self._make_group(inputs_tab, "Flowline Source: File or Source Flowlines")
        self._add_path(group, "flowline", "Flowline file", "", file_filter=vector_files)
        self._add_multi_files(group, "source_flowlines", "Source flowline files", _default("source_flowlines"), vector_files)

        group = self._make_group(run_tab, "Core")
        self._add_combo(group, "mapper", "Mapper", Mapper.list_names(), default=_default("mapper"))
        saved_streamflow_source = _streamflow_label(settings.value("streamflow_source", _default("streamflow_source")))
        self.streamflow_source = self._add_combo(group, "streamflow_source", "Streamflow source", STREAMFLOW_SOURCE_LABELS.keys(), default=saved_streamflow_source)
        self.nwm_api_key = self._add_line(group, "nwm_api_key", "NWM API key", "", placeholder="Required for NWM sources")
        self.nwm_api_key_label = group.itemAtPosition(group.rowCount() - 1, 0).widget()
        self._add_combo(group, "floodmap_mode", "Floodmap mode", [mode.value for mode in FloodMapMode], default=_default("floodmap_mode"))
        self._add_multi_files(group, "user_flow_files", "User flow files", _default("user_flow_files"), "CSV Files (*.csv);;All Files (*)")
        self._add_return_periods(group, "return_periods", "Return periods")
        self._add_path(group, "reanalysis_file", "Reanalysis file", "", file_filter="CSV Files (*.csv);;All Files (*)")

        group = self._make_group(run_tab, "Forecast")
        self._add_line(group, "forensic_forecast_date", "Forensic forecast date", "", placeholder="YYYYMMDD")
        self._add_combo(group, "forensic_forecast_hour", "Forensic forecast hour", [""] + [f"{hour:02d}" for hour in range(24)], default="")
        self._add_spin(group, "age_of_forecast_days", "Forecast age days", _default("age_of_forecast_days"), 1, 365)
        self._add_combo(group, "geoglows_vpu", "GEOGLOWS VPU", ["", "704", "702", "703", "715", "714", "706", "713", "712", "709"], default="")

        group = self._make_group(run_tab, "Execution")
        self._add_checkbox(group, "parallel", "Run workspaces in parallel", _default("parallel"))
        self._add_checkbox(group, "profile", "Print profiling stats", _default("profile"))
        self._add_checkbox(group, "quiet", "Reduce model output", _default("quiet"))
        self._add_checkbox(group, "use_yaml", "Write mapper configs as YAML", _default("use_yaml"))

        group = self._make_group(hydro_tab, "DEM and Stream Network")
        self._add_checkbox(group, "clean_dem", "Clean DEM", _default("clean_dem"))
        self._add_checkbox(group, "overwrite", "Overwrite existing products", _default("overwrite"))
        self._add_checkbox(group, "move_stream_network_to_thalweg", "Move stream network to thalweg", _default("move_stream_network_to_thalweg"))
        self._add_checkbox(group, "burn_streams", "Burn streams into DEM", _default("burn_streams"))
        self._add_checkbox(group, "project_to_utm", "Project everything to UTM", _default("project_to_utm"))
        self._add_checkbox(group, "streams_as_parquet", "Write streams as parquet", _default("streams_as_parquet"))
        self._add_checkbox(group, "use_parquet", "Use parquet intermediates", _default("use_parquet"))
        self._add_line(group, "new_strm_threshold_km2", "New stream threshold km2", _default("new_strm_threshold_km2"), parser="float")
        self._add_line(group, "min_match_score", "Minimum match score", "", parser="float")
        self._add_line(group, "q_baseflow_threshold", "Baseflow threshold", "", parser="float")

        group = self._make_group(hydro_tab, "Stream and Lake Filters")
        self._add_line(group, "StrmOrder_Field", "Stream order field", "")
        self._add_line(group, "StrmOrder_Lower", "Stream order lower", "", parser="int")
        self._add_line(group, "StrmOrder_Upper", "Stream order upper", "", parser="int")
        self._add_path(group, "lake_filter_json", "Lake filter JSON", "", file_filter="JSON Files (*.json);;All Files (*)")
        self._add_path(group, "lakes", "Lakes file", "", file_filter=vector_files)

        group = self._make_group(hydro_tab, "Land Cover")
        self._add_multi_files(group, "land_cover_cache", "Land cover cache files", _default("land_cover_cache"), "Raster Files (*.tif *.tiff *.vrt);;All Files (*)")
        self._add_line(group, "land_watervalue", "Land water value", _default("land_watervalue"), parser="int")
        self._add_checkbox(group, "flood_waterlc_and_strm_cells", "Flood water land-cover and stream cells", _default("flood_waterlc_and_strm_cells"))

        group = self._make_group(bathy_tab, "Bathymetry")
        self._add_checkbox(group, "disable_bathymetry", "Disable bathymetry", _default("disable_bathymetry"))
        self._add_checkbox(group, "bathy_use_banks", "Use bathymetry banks", _default("bathy_use_banks"))
        self._add_checkbox(group, "find_banks_based_on_landcover", "Find banks from land cover", _default("find_banks_based_on_landcover"))
        self._add_checkbox(group, "use_specified_depth_for_bathy_mask", "Use specified depths for bathymetry mask", _default("use_specified_depth_for_bathy_mask"))
        self._add_line(group, "specify_depths_for_bathy_mask", "Bathymetry mask depths", ", ".join(str(value) for value in _default("specify_depths_for_bathy_mask")), parser="float_list", placeholder="0.1 or 0.1, 0.3")
        self._add_combo(group, "specified_bathyflow_field", "Bathymetry flow field", FLOW_FIELD_OPTIONS, default=_default("specified_bathyflow_field"), editable=True)
        self._add_combo(group, "specified_highflow_field", "High flow field", FLOW_FIELD_OPTIONS, default=_default("specified_highflow_field"), editable=True)
        self._add_checkbox(group, "create_reach_average_curve_file", "Create reach-average curve file", _default("create_reach_average_curve_file"))
        self._add_checkbox(group, "use_power_laws_for_bathymetry", "Use power laws for bathymetry", _default("use_power_laws_for_bathymetry"))
        self._add_line(group, "area_m2_field", "Area m2 field", _default("area_m2_field"))
        self._add_line(group, "area_km2_field", "Area km2 field", _default("area_km2_field"))
        self._add_path(group, "mannings_text_file", "Manning's n text file", "", file_filter="Text Files (*.txt);;All Files (*)")

        group = self._make_group(bathy_tab, "Bathymetry Arguments")
        self._add_dict(group, "bathy_args", "Arguments", DEFAULT_BATHY_ARGS)

        group = self._make_group(outputs_tab, "Products")
        self._add_checkbox(group, "make_fist_inputs", "Make FIST inputs", _default("make_fist_inputs"))
        self._add_checkbox(group, "make_vdt", "Make VDT database", _default("make_vdt"))
        self._add_checkbox(group, "make_curvefile", "Make curve file", _default("make_curvefile"))
        self._add_checkbox(group, "make_ap_database", "Make AP database", _default("make_ap_database"))
        self._add_checkbox(group, "make_cross_section_file", "Make cross-section file", _default("make_cross_section_file"))
        self._add_checkbox(group, "make_depth_maps", "Make depth maps", _default("make_depth_maps"))
        self._add_checkbox(group, "make_velocity_maps", "Make velocity maps", _default("make_velocity_maps"))
        self._add_checkbox(group, "make_wse_maps", "Make WSE maps", _default("make_wse_maps"))
        self._add_checkbox(group, "estimate_consequences", "Estimate consequences", _default("estimate_consequences"))

        group = self._make_group(outputs_tab, "Flood Maps")
        self._add_checkbox(group, "overwrite_floodmaps", "Overwrite flood maps", _default("overwrite_floodmaps"))
        self._add_checkbox(group, "remove_old_forecast_files", "Remove old forecast files", _default("remove_old_forecast_files"))
        self._add_line(group, "floodmap_identifier", "Flood map identifier", "")
        self._add_combo(group, "vdt_file_extension", "VDT file extension", ["txt", "csv", "parquet"], default=_default("vdt_file_extension"))

        group = self._make_group(outputs_tab, "Floodmap Arguments")
        self._add_dict(group, "floodmap_args", "Arguments", DEFAULT_FLOODMAP_ARGS)

        group = self._make_group(advanced_tab, "FLDPLN")
        self._add_line(group, "fldpln_dh", "Depth interval", _default("fldpln_dh"), parser="float")
        self._add_line(group, "fldpln_min_depth", "Minimum depth", _default("fldpln_min_depth"), parser="float")
        self._add_line(group, "fldpln_max_depth", "Maximum depth", _default("fldpln_max_depth"), parser="float")
        self._add_line(group, "fldpln_max_wse_rise", "Maximum WSE rise", _default("fldpln_max_wse_rise"), parser="float")
        self._add_checkbox(group, "fldpln_keep_spilling", "Keep spilling cells", _default("fldpln_keep_spilling"))
        self._add_checkbox(group, "fldpln_parallel", "Run FLDPLN in parallel", _default("fldpln_parallel"))

        group = self._make_group(advanced_tab, "Runtime and Output")
        self._add_combo(group, "compression", "Raster compression", ["LZW", "DEFLATE", "ZSTD", "NONE"], default=_default("compression"), editable=True)
        self._add_checkbox(group, "raise_errors_if_nothing_in_domain", "Raise if nothing is in domain", _default("raise_errors_if_nothing_in_domain"))
        self._add_line(group, "exclude", "Exclude stream IDs", "", parser="int_list", placeholder="123, 456")

    def _wire_map_refresh(self):
        map_keys = {
            "dem",
            "dem_dir",
            "flowline",
            "source_dems",
            "land_cover_cache",
            "user_flow_files",
        }
        for key in map_keys:
            widget = self.input_fields.get(key)
            if isinstance(widget, PathPicker):
                widget.line_edit.textChanged.connect(self.schedule_input_layer_refresh)
            elif isinstance(widget, MultiFilePicker):
                widget.changed_signal.connect(self.schedule_input_layer_refresh)

        self.schedule_input_layer_refresh()

    def schedule_input_layer_refresh(self):
        if hasattr(self, "map_refresh_timer"):
            self.map_refresh_timer.start(500)

    def refresh_input_layers(self):
        specs = []
        raster_filter = {".tif", ".tiff", ".vrt", ".img"}
        vector_filter = {".gpkg", ".shp", ".gdb", ".parquet", ".geoparquet", ".geojson", ".json"}

        dem_dir = self.input_fields["dem_dir"].text().strip()
        dem_file = self.input_fields["dem"].text().strip()
        if dem_file:
            specs.append(("DEM file", dem_file))
        if dem_dir:
            for path in MapLayerLoader.discover_rasters(dem_dir)[:50]:
                specs.append((f"DEM directory/{path.name}", path))

        for path in self.input_fields["source_dems"].values():
            specs.append((f"Source DEM/{Path(path).name}", path))
        flowline = self.input_fields["flowline"].text().strip()
        if flowline:
            specs.append((f"Flowline/{Path(flowline).name}", flowline))
        for path in self.input_fields["land_cover_cache"].values():
            specs.append((f"Land cover/{Path(path).name}", path))
        for path in self.input_fields["user_flow_files"].values():
            specs.append((f"User flow file/{Path(path).name}", path))

        layers = []
        for index, (name, path) in enumerate(specs):
            suffix = Path(path).suffix.lower()
            if suffix in raster_filter | vector_filter or Path(path).exists():
                layers.append(MapLayerLoader.make_layer(name, path, "Inputs", index))
        self.set_map_group_layers("Inputs", layers)

    def load_output_layers(self):
        try:
            params = self._get_params()
        except Exception as exc:
            self.log_message(f"[WARNING] Could not collect output path settings: {exc}")
            return

        output_dir = params.get("output_dir")
        watershed_name = params.get("name")
        if not output_dir or not watershed_name:
            self.log_message("[WARNING] Output directory and watershed name are required to load outputs.")
            return

        paths = MapLayerLoader.discover_outputs(output_dir, watershed_name)
        stream_layers = self.make_generated_stream_layers(params)
        layers = stream_layers + [
            MapLayerLoader.make_layer(str(path.relative_to(Path(output_dir) / watershed_name)), path, "Outputs", index)
            for index, path in enumerate(paths[:150])
            if path not in {layer["path"] for layer in stream_layers}
        ]
        self.set_map_group_layers("Outputs", layers)
        self.side_tabs.setCurrentIndex(1)
        self.log_message(f"[INFO] Loaded {len(layers)} output layers.")

    def make_generated_stream_layers(self, params):
        try:
            from nencarta.workspace import Workspace

            configs = NencartaConfig(dict(params))
            dem_inputs = []
            if configs.dem_dir:
                dem_inputs = list(Path(configs.dem_dir).glob(configs.dem_filter))
            elif configs.dem:
                dem_inputs = [Path(configs.dem)]
            elif configs.source_dems and configs.bbox:
                dem_inputs = [None]

            layers = []
            for index, dem in enumerate(dem_inputs[:50]):
                workspace = Workspace(configs, dem)
                if configs.move_stream_network_to_thalweg and workspace.new_StrmShp_matched.exists():
                    path = workspace.new_StrmShp_matched
                    name = f"Generated streams matched/{workspace.FileName}"
                elif workspace.DEM_StrmShp.exists():
                    path = workspace.DEM_StrmShp
                    name = f"Generated streams/{workspace.FileName}"
                else:
                    continue
                layers.append(MapLayerLoader.make_layer(name, path, "Outputs", index))
            return layers
        except Exception as exc:
            self.log_message(f"[WARNING] Could not inspect generated stream outputs: {exc}")
            return []

    def set_map_group_layers(self, group, layers):
        self.layer_groups[group] = layers
        self.render_map_layers()

    def clear_map_layers(self):
        self.layer_groups = {"Inputs": [], "Outputs": []}
        self.layer_visibility.clear()
        self.stop_basemap_workers()
        self.render_map_layers()

    def render_map_layers(self):
        for item in self.map_scene.items():
            self.map_scene.removeItem(item)
        self.clear_layer_controls()
        has_layers = False
        layer_rect = None

        for group_name, layers in self.layer_groups.items():
            if not layers:
                continue
            has_layers = True
            group_label = QLabel(group_name)
            group_label.setObjectName("LayerGroupLabel")
            self.layer_controls_layout.insertWidget(self.layer_controls_layout.count() - 1, group_label)
            for layer in layers:
                key = f"{group_name}:{layer['path']}:{layer['name']}"
                visible = self.layer_visibility.get(key, True)
                checkbox = QCheckBox(layer["name"])
                checkbox.setChecked(visible)
                checkbox.setToolTip(str(layer["path"]))
                checkbox.toggled.connect(lambda checked, item=layer["item"], layer_key=key: self.set_layer_visible(layer_key, item, checked))

                stats_button = QPushButton("Stats")
                stats_button.setCheckable(True)
                stats_button.setMaximumWidth(70)
                header = QHBoxLayout()
                header.addWidget(checkbox, 1)
                header.addWidget(stats_button)

                stats_text = QPlainTextEdit(layer["stats"])
                stats_text.setReadOnly(True)
                stats_text.setVisible(False)
                stats_text.setMaximumHeight(150)
                stats_button.toggled.connect(stats_text.setVisible)

                layer_widget = QWidget()
                layer_layout = QVBoxLayout(layer_widget)
                layer_layout.setContentsMargins(0, 0, 0, 0)
                layer_layout.setSpacing(4)
                layer_layout.addLayout(header)
                layer_layout.addWidget(stats_text)
                self.layer_controls_layout.insertWidget(self.layer_controls_layout.count() - 1, layer_widget)

                if layer["item"] is not None:
                    layer["item"].setVisible(visible)
                    self.map_scene.addItem(layer["item"])
                    if visible:
                        item_rect = layer["item"].sceneBoundingRect()
                        layer_rect = item_rect if layer_rect is None else layer_rect.united(item_rect)

        if not has_layers:
            empty_label = QLabel("No map layers loaded.")
            empty_label.setObjectName("LayerEmptyLabel")
            self.layer_controls_layout.insertWidget(self.layer_controls_layout.count() - 1, empty_label)

        self.current_layer_rect = layer_rect
        self.render_basemap(layer_rect)
        if self.map_scene.items():
            self.zoom_to_map_layers(refresh_basemap=False)

    def render_basemap(self, layer_rect):
        self.remove_basemap_items()
        if not self.google_hybrid_basemap.isChecked() or layer_rect is None or layer_rect.isNull():
            return

        self.basemap_request_id += 1
        request_id = self.basemap_request_id
        bbox = self.scene_rect_to_bbox(layer_rect)
        tile_specs = MapLayerLoader.google_hybrid_tile_specs(bbox)
        if not tile_specs:
            self.log_basemap_warning_once()
            return

        missing_tiles = []
        for zoom, x, y in tile_specs:
            key = (zoom, x, y)
            data = self.basemap_tile_cache.get(key)
            if data is None:
                missing_tiles.append((zoom, x, y))
            else:
                self.add_basemap_tile(request_id, zoom, x, y, data)

        if missing_tiles:
            worker = BasemapTileWorker(request_id, missing_tiles)
            worker.tile_signal.connect(self.add_basemap_tile)
            worker.finished_signal.connect(self.basemap_worker_finished)
            worker.finished.connect(lambda worker=worker: self.cleanup_basemap_worker(worker))
            self.basemap_workers.append(worker)
            worker.start()

    def refresh_basemap_for_view(self):
        if not self.google_hybrid_basemap.isChecked() or self.current_layer_rect is None:
            return
        visible_rect = self.map_view.mapToScene(self.map_view.viewport().rect()).boundingRect()
        if not visible_rect.isValid() or visible_rect.isNull():
            visible_rect = self.current_layer_rect
        else:
            visible_rect = visible_rect.intersected(self.current_layer_rect)
            if visible_rect.isNull():
                visible_rect = self.current_layer_rect
        self.render_basemap(visible_rect)

    def schedule_basemap_refresh(self):
        if hasattr(self, "basemap_refresh_timer") and self.google_hybrid_basemap.isChecked():
            self.basemap_refresh_timer.start(850)

    def remove_basemap_items(self):
        self.stop_basemap_workers()
        for item in self.basemap_items:
            if item.scene() is self.map_scene:
                self.map_scene.removeItem(item)
        self.basemap_items = []

    def add_basemap_tile(self, request_id, zoom, x, y, data):
        if request_id != self.basemap_request_id:
            return
        key = (zoom, x, y)
        self.basemap_tile_cache[key] = data
        while len(self.basemap_tile_cache) > 2048:
            self.basemap_tile_cache.pop(next(iter(self.basemap_tile_cache)))

        item = MapLayerLoader.make_google_tile_item_from_data(x, y, zoom, data)
        if item is None:
            return
        self.basemap_items.append(item)
        self.map_scene.addItem(item)

    def basemap_worker_finished(self, request_id, loaded_count):
        if request_id == self.basemap_request_id and loaded_count == 0 and not self.basemap_items:
            self.log_basemap_warning_once()

    def cleanup_basemap_worker(self, worker):
        if worker in self.basemap_workers:
            self.basemap_workers.remove(worker)
        worker.deleteLater()

    def stop_basemap_workers(self):
        self.basemap_request_id += 1
        for worker in list(self.basemap_workers):
            if worker.isRunning():
                worker.requestInterruption()

    def log_basemap_warning_once(self):
        if not self._basemap_warning_shown:
            self.log_message("[WARNING] Google Hybrid basemap tiles could not be loaded.")
            self._basemap_warning_shown = True

    @staticmethod
    def scene_rect_to_bbox(rect):
        return (
            rect.left(),
            -rect.bottom(),
            rect.right(),
            -rect.top(),
        )

    def clear_layer_controls(self):
        while self.layer_controls_layout.count() > 1:
            item = self.layer_controls_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()

    def set_layer_visible(self, key, item, visible):
        self.layer_visibility[key] = visible
        if item is not None:
            item.setVisible(visible)

    def zoom_to_map_layers(self, refresh_basemap=True):
        rect = self.map_scene.itemsBoundingRect()
        if rect.isValid() and not rect.isNull():
            margin_x = max(rect.width() * 0.08, 0.01)
            margin_y = max(rect.height() * 0.08, 0.01)
            self.map_view.fitInView(rect.adjusted(-margin_x, -margin_y, margin_x, margin_y), Qt.KeepAspectRatio)
            if refresh_basemap:
                self.schedule_basemap_refresh()

    def _wire_visibility(self):
        def toggle_nwm_api_key(value):
            source = STREAMFLOW_SOURCE_LABELS.get(value, value)
            is_nwm = str(source).upper().startswith("NWM")
            self.nwm_api_key.setVisible(is_nwm)
            self.nwm_api_key_label.setVisible(is_nwm)

        self.streamflow_source.currentTextChanged.connect(toggle_nwm_api_key)
        toggle_nwm_api_key(self.streamflow_source.currentText())

    def _parse_widget_value(self, key, widget):
        parser = self.field_parsers.get(key, "text")
        if isinstance(widget, PathPicker):
            value = widget.text().strip()
        elif isinstance(widget, QLineEdit):
            value = widget.text().strip()
        elif isinstance(widget, QPlainTextEdit):
            value = widget.toPlainText().strip()
        elif isinstance(widget, QCheckBox):
            return widget.isChecked()
        elif isinstance(widget, QComboBox):
            value = widget.currentText().strip()
        elif isinstance(widget, QSpinBox):
            return widget.value()
        elif isinstance(widget, DictTable):
            return widget.to_dict()
        elif isinstance(widget, MultiFilePicker):
            return widget.values()
        elif isinstance(widget, ReturnPeriodPicker):
            return widget.values()
        else:
            return None

        if parser in {"text", "path"}:
            return value or None
        if parser == "str_list":
            return _split_list(value)
        if parser == "float_list":
            return _parse_float_list(value)
        if parser == "int_list":
            return _parse_int_list(value)
        if parser == "bbox":
            return _parse_bbox(value)
        if value == "":
            return None
        if parser == "int":
            return int(float(value))
        if parser == "float":
            return float(value)
        if parser == "bool":
            return bool(value)
        return value

    def _get_params(self):
        params = {}
        for key, widget in self.input_fields.items():
            params[key] = self._parse_widget_value(key, widget)

        params["name"] = params.pop("watershed_name")
        streamflow_label = params.get("streamflow_source")
        if streamflow_label:
            params["streamflow_source"] = STREAMFLOW_SOURCE_LABELS.get(streamflow_label, streamflow_label)

        clean_params = {}
        for key, value in params.items():
            if value is None or value == "" or value == []:
                continue
            clean_params[key] = value
        return clean_params

    def _validate_params(self, params):
        if not params.get("name"):
            raise ValueError("Watershed name is required.")
        if not params.get("output_dir"):
            raise ValueError("Output directory is required.")
        if not any(params.get(key) for key in ("flowline", "source_flowlines")):
            raise ValueError("Provide a flowline file or source flowline files.")
        if not any(params.get(key) for key in ("dem", "dem_dir", "source_dems")):
            raise ValueError("Provide a DEM file, DEM directory, or source DEM files.")
        if params.get("streamflow_source", "").upper().startswith("NWM") and not params.get("nwm_api_key"):
            raise ValueError("NWM API key is required for NWM streamflow sources.")
        if params.get("floodmap_mode") == FloodMapMode.USER.value and not params.get("user_flow_files"):
            raise ValueError("User flow files are required when floodmap mode is user.")
        if params.get("floodmap_mode") == FloodMapMode.RETURN_PERIOD.value and not params.get("return_periods"):
            raise ValueError("Return periods are required when floodmap mode is return_period.")

        NencartaConfig(dict(params))

    def start_simulation(self):
        if self.worker_thread and self.worker_thread.isRunning():
            QMessageBox.warning(self, "Running", "A simulation is already running.")
            return

        self.log_text.clear()
        self.log_text.append("[INFO] Collecting parameters...")
        self.run_button.setEnabled(False)
        self.run_button.setText("Running...")

        try:
            params = self._get_params()
            self._validate_params(params)
            if params.get("dem") and params.get("dem_dir"):
                self.log_text.append(
                    "[WARNING] Both DEM directory and DEM file are set. "
                    "process_watershed uses the DEM directory first."
                )
            self.save_settings(params)
        except Exception as e:
            self.show_error(f"Parameter error: {e}")
            return

        self.worker_thread = WorkerThread(params)
        self.worker_thread.finished_signal.connect(self.display_results)
        self.worker_thread.log_signal.connect(self.log_message)
        self.worker_thread.error_signal.connect(self.show_error)
        self.worker_thread.start()

    def log_message(self, message):
        message = str(message)
        if message.startswith("\r"):
            if self._last_log_message_was_progress:
                cursor = self.log_text.textCursor()
                cursor.movePosition(QTextCursor.End)
                cursor.select(QTextCursor.LineUnderCursor)
                cursor.removeSelectedText()
                cursor.insertText(message[1:])
            else:
                self.log_text.append(message[1:])
            self._last_log_message_was_progress = True
        else:
            self.log_text.append(message)
            self._last_log_message_was_progress = False
        self.log_text.verticalScrollBar().setValue(self.log_text.verticalScrollBar().maximum())

    def display_results(self, watershed_name: str):
        self.log_message("[COMPLETED] Simulation finished.")
        self.load_output_layers()
        self.run_button.setEnabled(True)
        self.run_button.setText("Start Simulation")

    def show_error(self, message):
        summary = str(message)
        details = None
        marker = "\n\nDETAILS:\n"
        if marker in summary:
            summary, details = summary.split(marker, 1)

        self.log_text.append(f"[ERROR] {summary}")

        dlg = QMessageBox(self)
        dlg.setIcon(QMessageBox.Critical)
        dlg.setWindowTitle("Simulation Error")
        dlg.setText("An error occurred.")
        dlg.setInformativeText(summary)
        if details:
            dlg.setDetailedText(details)
        dlg.exec_()

        self.run_button.setEnabled(True)
        self.run_button.setText("Start Simulation")

    def closeEvent(self, event):
        self.stop_basemap_workers()
        for worker in list(self.basemap_workers):
            if worker.isRunning():
                worker.wait(1200)
        super().closeEvent(event)

    def save_settings(self, params=None):
        if params is None:
            for key, widget in self.input_fields.items():
                setting_key = "watershed_name" if key == "watershed_name" else key
                if isinstance(widget, PathPicker):
                    settings.setValue(setting_key, widget.text().strip())
                elif isinstance(widget, QLineEdit):
                    settings.setValue(setting_key, widget.text().strip())
                elif isinstance(widget, QPlainTextEdit):
                    settings.setValue(setting_key, widget.toPlainText().strip())
                elif isinstance(widget, QCheckBox):
                    settings.setValue(setting_key, widget.isChecked())
                elif isinstance(widget, QComboBox):
                    settings.setValue(setting_key, widget.currentText().strip())
                elif isinstance(widget, QSpinBox):
                    settings.setValue(setting_key, widget.value())
                elif isinstance(widget, DictTable):
                    settings.setValue(setting_key, json.dumps(widget.to_dict()))
                elif isinstance(widget, MultiFilePicker):
                    settings.setValue(setting_key, json.dumps(widget.values()))
                elif isinstance(widget, ReturnPeriodPicker):
                    settings.setValue(setting_key, json.dumps(widget.values()))
            return

        for key, value in params.items():
            if key == "name":
                settings.setValue("watershed_name", value)
            elif key == "streamflow_source":
                settings.setValue(key, _streamflow_label(value))
            elif isinstance(value, (dict, list)):
                settings.setValue(key, json.dumps(value))
            else:
                settings.setValue(key, value)


def apply_theme(app):
    app.setStyle("Fusion")
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor("#f6f7f9"))
    palette.setColor(QPalette.WindowText, QColor("#1f2933"))
    palette.setColor(QPalette.Base, QColor("#ffffff"))
    palette.setColor(QPalette.AlternateBase, QColor("#f3f5f7"))
    palette.setColor(QPalette.ToolTipBase, QColor("#ffffff"))
    palette.setColor(QPalette.ToolTipText, QColor("#1f2933"))
    palette.setColor(QPalette.Text, QColor("#1f2933"))
    palette.setColor(QPalette.Button, QColor("#ffffff"))
    palette.setColor(QPalette.ButtonText, QColor("#1f2933"))
    palette.setColor(QPalette.BrightText, QColor("#b42318"))
    palette.setColor(QPalette.Link, QColor("#2563eb"))
    palette.setColor(QPalette.Highlight, QColor("#2563eb"))
    palette.setColor(QPalette.HighlightedText, QColor("#ffffff"))
    app.setPalette(palette)

    app.setStyleSheet(
        """
        QWidget {
            font-family: Segoe UI, Arial, sans-serif;
            font-size: 10pt;
            color: #1f2933;
        }
        QMainWindow, QWidget#centralWidget {
            background: #f6f7f9;
        }
        QLabel#AppTitle {
            font-size: 22pt;
            font-weight: 700;
            color: #111827;
        }
        QLabel#AppSubtitle {
            color: #5f6b7a;
            font-size: 10pt;
        }
        QGroupBox {
            background: #ffffff;
            border: 1px solid #d9dee7;
            border-radius: 8px;
            margin-top: 12px;
            font-weight: 600;
        }
        QGroupBox::title {
            subcontrol-origin: margin;
            left: 12px;
            padding: 0 4px;
        }
        QLabel#LayerGroupLabel {
            color: #111827;
            font-weight: 700;
            padding-top: 6px;
        }
        QGraphicsView {
            background: #eef2f7;
            border: 1px solid #cfd6df;
            border-radius: 6px;
        }
        QListWidget {
            background: #ffffff;
            border: 1px solid #cfd6df;
            border-radius: 6px;
            padding: 4px;
        }
        QLineEdit, QPlainTextEdit, QTextEdit, QComboBox, QSpinBox, QTableWidget {
            background: #ffffff;
            border: 1px solid #cfd6df;
            border-radius: 6px;
            padding: 6px;
            selection-background-color: #2563eb;
        }
        QLineEdit:focus, QPlainTextEdit:focus, QTextEdit:focus, QComboBox:focus, QSpinBox:focus {
            border: 1px solid #2563eb;
        }
        QTextEdit {
            background: #111827;
            color: #e5e7eb;
            border: 1px solid #1f2937;
        }
        QPushButton {
            background: #ffffff;
            border: 1px solid #cfd6df;
            border-radius: 6px;
            padding: 7px 12px;
        }
        QPushButton:hover {
            background: #f3f5f7;
        }
        QPushButton#PrimaryButton {
            background: #2563eb;
            color: #ffffff;
            border: 1px solid #2563eb;
            font-weight: 600;
            padding: 8px 18px;
        }
        QPushButton#PrimaryButton:hover {
            background: #1d4ed8;
        }
        QPushButton:disabled {
            color: #8a95a3;
            background: #e7ebf0;
            border-color: #d7dde5;
        }
        QTabWidget::pane {
            border: 1px solid #d9dee7;
            border-radius: 8px;
            background: #ffffff;
        }
        QTabBar::tab {
            background: #e7ebf0;
            border: 1px solid #d9dee7;
            border-bottom: none;
            border-top-left-radius: 6px;
            border-top-right-radius: 6px;
            padding: 8px 12px;
            margin-right: 2px;
        }
        QTabBar::tab:selected {
            background: #ffffff;
            color: #111827;
        }
        QScrollArea {
            border: none;
            background: transparent;
        }
        QCheckBox {
            spacing: 8px;
        }
        QTableWidget {
            gridline-color: #d9dee7;
        }
        QHeaderView::section {
            background: #f3f5f7;
            border: 0;
            border-bottom: 1px solid #d9dee7;
            padding: 6px;
            font-weight: 600;
        }
        """
    )


def run_gui():
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)

    apply_theme(app)
    window = FloodSimulationGUI()
    app.aboutToQuit.connect(window.save_settings)
    window.show()
    sys.exit(app.exec_())
