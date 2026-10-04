"""
Burning the source streams into the DEM, and moving them onto its D8 flow paths (nencarta.tasks.move_streams).

The burn has one job the D8 flow paths depend on: every burned river has to keep draining downstream. Where it
doesn't, fill turns the burned channel above the blockage into a flat pool and the flow paths across it run in
straight lines -- on the N38W085 tile that is what flattened 30 km of the Ohio.
"""
import heapq

import numpy as np
import pandas as pd
import geopandas as gpd
import pytest
from osgeo import gdal, osr
from shapely.geometry import LineString, box

from nencarta.core.configs import NencartaConfig
from nencarta.core.raster import Raster
from nencarta.workspace import Workspace
from nencarta.tasks.move_streams import (
    _drop_reaches,
    _get_stream_raster,
    _reaches_stuck_in_lakes,
    burn_streams_and_move_streams,
    burn_streams_into_dem,
)
from nencarta.tasks.stream_conflation import D8_OFFSETS

CELL = 30.0
X0, Y0 = 500_000.0, 4_500_000.0
GT = (X0, CELL, 0.0, Y0, 0.0, -CELL)


def utm_wkt() -> str:
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32615)
    return srs.ExportToWkt()


def centre(r, c):
    return X0 + (c + 0.5) * CELL, Y0 - (r + 0.5) * CELL


def mem_dataset(array: np.ndarray, nodata: float = -9999.0) -> gdal.Dataset:
    ds = gdal.GetDriverByName("MEM").Create("", array.shape[1], array.shape[0], 1, gdal.GDT_Float32)
    ds.SetGeoTransform(GT)
    ds.SetProjection(utm_wkt())
    ds.GetRasterBand(1).WriteArray(array)
    ds.GetRasterBand(1).SetNoDataValue(nodata)
    return ds


def write_raster(path, array, dtype=gdal.GDT_Float32, nodata=None):
    ds = gdal.GetDriverByName("GTiff").Create(str(path), array.shape[1], array.shape[0], 1, dtype)
    ds.SetGeoTransform(GT)
    ds.SetProjection(utm_wkt())
    ds.GetRasterBand(1).WriteArray(array)
    if nodata is not None:
        ds.GetRasterBand(1).SetNoDataValue(nodata)
    ds = None


def filled(dem: np.ndarray) -> np.ndarray:
    """A priority-flood depression fill draining off every edge (what WhiteboxTools does before D8)."""
    nrows, ncols = dem.shape
    out = dem.astype(np.float64).copy()
    done = np.zeros(dem.shape, bool)
    heap = []
    for r in range(nrows):
        for c in range(ncols):
            if r in (0, nrows - 1) or c in (0, ncols - 1):
                heapq.heappush(heap, (out[r, c], r, c))
                done[r, c] = True
    while heap:
        z, r, c = heapq.heappop(heap)
        for dr, dc in D8_OFFSETS.values():
            nr, nc = r + dr, c + dc
            if 0 <= nr < nrows and 0 <= nc < ncols and not done[nr, nc]:
                done[nr, nc] = True
                out[nr, nc] = max(out[nr, nc], z)
                heapq.heappush(heap, (out[nr, nc], nr, nc))
    return out


# ---------------------------------------------------------------------------------------------------------------
# The burn
# ---------------------------------------------------------------------------------------------------------------

@pytest.fixture
def flat_river():
    """
    A 30 x 40 tile with a river's flattened water surface at 141.2 m across rows 10 to 20 and land at 150 m. The
    river (row 15, flowing east off the edge) is split at column 20, where a tributary from the north joins it,
    and the network lists the tributary last, so drawn in table order it would own the cell where the three meet.
    A bank cell at 135 m near the river's upstream end drags the whole burned channel below the water surface,
    as FABDEM's flattened rivers do.
    """
    dem = np.full((30, 40), 150.0, np.float32)
    dem[10:21, :] = 141.2
    dem[14, 1] = 135.0
    for r in range(0, 10):
        dem[r, 18:23] = 146.0 - 0.4 * r  # the tributary's valley
    upper = LineString([centre(15, c) for c in range(0, 21)][::-1])   # drawn downstream to upstream, as TDX-Hydro
    lower = LineString([centre(15, c) for c in range(20, 40)][::-1])
    tributary = LineString([centre(r, 20) for r in range(0, 16)][::-1])
    gdf = gpd.GeoDataFrame({"LINKNO": [3, 1, 2], "DSLINKNO": [-1, 3, 3], "strmOrder": [2, 2, 1],
                            "DSContArea": [900.0, 850.0, 20.0]},
                           geometry=[lower, upper, tributary], crs="EPSG:32615")
    channel_mask = np.zeros(dem.shape, bool)
    channel_mask[10:21, :] = True
    return dem, gdf, channel_mask


def burn(dem, gdf, channel_mask, order_col="strmOrder", area_col="DSContArea"):
    dem = dem.copy()
    dem_ds = mem_dataset(dem)
    streams = _get_stream_raster(dem_ds, gdf, "LINKNO", order_col, area_col)
    burn_streams_into_dem(channel_mask.copy(), dem, streams, dem_ds, gdf, None, "LINKNO", "DSLINKNO", -9999.0)
    return dem, streams


def test_the_larger_reach_owns_the_cell_where_reaches_meet(flat_river):
    dem, gdf, _ = flat_river
    streams = _get_stream_raster(mem_dataset(dem), gdf, "LINKNO", "strmOrder", "DSContArea")
    assert streams[15, 20] == 3  # the river below the confluence, not the tributary listed after it

    # Without order or area to go on, the table order still decides
    streams = _get_stream_raster(mem_dataset(dem), gdf, "LINKNO")
    assert streams[15, 20] == 2


def test_a_burned_river_keeps_draining_through_a_confluence(flat_river):
    dem, gdf, channel_mask = flat_river
    burned, streams = burn(dem, gdf, channel_mask)

    river = [(15, c) for c in range(0, 40)]
    profile = np.array([burned[cell] for cell in river])
    assert (np.diff(profile[1:]) <= 0).all(), profile
    # Nothing along the river is a pit that filling would flatten
    assert np.allclose(filled(burned)[15, :], burned[15, :])


def test_a_one_cell_gap_in_the_channel_mask_is_burned_through(flat_river):
    # The channel mask is made from the cleaned stream raster, which can leave out a cell the burn's own stream
    # raster has. The burn used to skip such a cell, leaving it at the water's surface: a dam.
    dem, gdf, channel_mask = flat_river
    channel_mask[15, 30] = False
    burned, _ = burn(dem, gdf, channel_mask)

    assert burned[15, 30] <= burned[15, 29]
    assert np.allclose(filled(burned)[15, :], burned[15, :])


def test_a_reach_with_nothing_up_or_downstream_is_burned():
    dem = np.full((20, 20), 150.0, np.float32)
    dem[:, 10] = 140.3
    gdf = gpd.GeoDataFrame({"LINKNO": [7], "DSLINKNO": [-1]},
                           geometry=[LineString([centre(r, 10) for r in range(0, 20)])], crs="EPSG:32615")
    burned, _ = burn(dem, gdf, np.zeros(dem.shape, bool), order_col=None, area_col=None)

    assert (burned[1:19, 10] <= 140.0).all()


# ---------------------------------------------------------------------------------------------------------------
# Lakes, for the mappers other than FLDPLN, and dropping reaches
# ---------------------------------------------------------------------------------------------------------------

def network(rows):
    return gpd.GeoDataFrame({"LINKNO": [r[0] for r in rows], "DSLINKNO": [r[1] for r in rows]},
                            geometry=[LineString([(0, 0), (1, 1)])] * len(rows))


def test_lake_reaches_go_unless_they_carry_a_river_through_the_lake():
    #   1 -> 2 (in the lake) -> 3         2 carries the river through: kept
    #   4 -> 5 (in the lake), an outlet   5 leads nowhere outside the lake: dropped
    #   6 (in the lake) -> 7              6 has nothing above it outside the lake: dropped
    streams = network([(1, 2), (2, 3), (3, -1), (4, 5), (5, -1), (6, 7), (7, -1)])
    raster = np.array([[1, 2, 3, 4, 5, 6, 7]], np.int32)
    lakes = np.array([[False, True, False, False, True, True, False]])

    assert _reaches_stuck_in_lakes(streams, raster, lakes, "LINKNO", "DSLINKNO") == {5, 6}


def test_dropping_a_reach_clears_its_cells_and_frees_the_reaches_above_it():
    streams = network([(1, 2), (2, 3), (3, -1)])
    raster = np.array([[1, 1, 2, 2, 3]], np.int32)

    streams, raster = _drop_reaches(streams, raster, {2}, "LINKNO", "DSLINKNO")

    assert list(streams.LINKNO) == [1, 3]
    assert streams.set_index("LINKNO").at[1, "DSLINKNO"] == -1
    assert raster.tolist() == [[1, 1, 0, 0, 3]]


# ---------------------------------------------------------------------------------------------------------------
# The whole step, through WhiteboxTools
# ---------------------------------------------------------------------------------------------------------------

def synthetic_watershed(tmp, mapper="Curve2Flood-FLDPLNpy", with_lake=True):
    """
    An 80 x 80 tile of V-shaped valleys: a river flowing east along row 40 and off the east edge, and a tributary
    from the north joining it at column 40. A lake covers the river from column 60 to 68, and reach 14 crosses it,
    with river reaches above and below it. The source network draws its lines a cell off the valley floors, and
    from the outlet up, as TDX-Hydro does.
    """
    rows, cols = np.mgrid[0:80, 0:80]
    main = 10 + 0.2 * (79 - cols) + 2.0 * np.abs(rows - 40)
    trib = np.where(rows <= 40, 10 + 0.2 * 39 + 0.2 * (40 - rows) + 2.0 * np.abs(cols - 40), np.inf)
    dem = np.minimum(main, trib).astype(np.float32)
    dem_path = tmp / "dem.tif"
    write_raster(dem_path, dem, nodata=-9999.0)

    def line(cells):
        return LineString([centre(r + 0.4, c - 0.3) for r, c in cells][::-1])

    streams = gpd.GeoDataFrame(
        {"LINKNO": [11, 12, 13, 14, 15], "DSLINKNO": [13, 13, 14, 15, -1], "strmOrder": [2, 1, 2, 2, 2],
         "DSContArea": [2.0e6, 1.0e6, 3.5e6, 4.0e6, 4.2e6]},
        geometry=[line([(40, c) for c in range(0, 41)]), line([(r, 40) for r in range(0, 41)]),
                  line([(40, c) for c in range(40, 57)]), line([(40, c) for c in range(56, 73)]),
                  line([(40, c) for c in range(72, 80)])],
        crs="EPSG:32615")
    config = {
        "name": "synthetic", "dem_dir": str(tmp), "flowline": str(tmp / "streams.gpkg"), "output_dir": str(tmp / "out"),
        "mapper": mapper,
        "streamflow_source": "GEOGLOWS", "stream_id_field": "LINKNO", "downstream_id_field": "DSLINKNO",
        "StrmOrder_Field": "strmOrder", "area_m2_field": "DSContArea", "move_stream_network_to_thalweg": True,
        "burn_streams": True, "new_strm_threshold_km2": 0.05, "clean_dem": False, "make_fist_inputs": False,
        "use_parquet": True, "streams_as_parquet": False, "short_file_names": True, "compression": "DEFLATE",
    }
    if with_lake:
        lake_path = tmp / "lakes.gpkg"
        x0, y1 = centre(34, 60)
        x1, y0 = centre(46, 68)
        gpd.GeoDataFrame({"Lake_id": [1]}, geometry=[box(x0, y0, x1, y1)], crs="EPSG:32615").to_file(lake_path)
        config["lakes"] = str(lake_path)
    workspace = Workspace(NencartaConfig(config), dem_path)
    workspace.DEM_StrmShp.parent.mkdir(parents=True, exist_ok=True)
    streams.to_file(workspace.DEM_StrmShp)
    water = np.zeros(dem.shape, np.uint8)
    water[40, :] = 1
    water[:41, 40] = 1
    workspace.bathy_water_mask.parent.mkdir(parents=True, exist_ok=True)
    write_raster(workspace.bathy_water_mask, water, gdal.GDT_Byte)
    return workspace


def read_d8_paths(workspace):
    flowdir = Raster(workspace.flowdir).read_array()
    raster = Raster(workspace.new_stream_raster).read_array()
    streams = gpd.read_file(workspace.new_StrmShp_matched)
    return flowdir, raster, streams


def test_burn_and_move_on_a_synthetic_watershed(workspace_tmp):
    workspace = synthetic_watershed(workspace_tmp)
    burn_streams_and_move_streams(workspace)

    flowdir, raster, streams = read_d8_paths(workspace)
    lakes = Raster(workspace.lake_raster).read_array().astype(bool)

    # FLDPLN doesn't map lakes: nothing runs inside one, the river above the lake ends at its shore, and the river
    # below it starts at the other shore
    assert {11, 12, 13, 15} <= set(streams.LINKNO)
    assert lakes.any() and not (raster.astype(bool) & lakes).any()
    ds = streams.set_index("LINKNO")["DSLINKNO"].to_dict()
    assert ds[11] == 13 and ds[12] == 13
    lake_cols = np.nonzero(lakes[40])[0]
    assert np.nonzero(raster == 13)[1].max() < lake_cols.min()
    assert np.nonzero(raster == 15)[1].min() > lake_cols.max()

    # Every reach is one D8 path, and its cells are the raster's
    for _, row in streams.iterrows():
        cells = [(int((Y0 - y) // CELL), int((x - X0) // CELL)) for x, y in row.geometry.coords]
        for a, b in zip(cells[:-1], cells[1:]):
            d_row, d_col = D8_OFFSETS[int(flowdir[a])]
            assert (a[0] + d_row, a[1] + d_col) == b
        assert all(raster[cell] == row.LINKNO for cell in cells)
    # Both rivers run down their valley floors, not where the source drew them
    assert set(np.unique(np.nonzero(raster == 11)[0])) == {40}
    assert set(np.unique(np.nonzero(raster == 12)[1])) == {40}

    table = pd.read_parquet(workspace.stream_info_file)
    assert set(table.stream_id) == set(streams.LINKNO)
    assert "topological_order" in streams.columns


def test_other_mappers_keep_a_lake_reach_that_carries_the_river(workspace_tmp):
    workspace = synthetic_watershed(workspace_tmp, mapper="Curve2Flood-Kernel Weighted")
    burn_streams_and_move_streams(workspace)

    _, raster, streams = read_d8_paths(workspace)
    # 14 crosses the lake with river above and below it, so the old rule keeps it, whole
    assert set(streams.LINKNO) == {11, 12, 13, 14, 15}
    assert streams.set_index("LINKNO")["DSLINKNO"].to_dict() == {11: 13, 12: 13, 13: 14, 14: 15, 15: -1}
    lakes = Raster(workspace.lake_raster).read_array().astype(bool)
    assert (raster[lakes] == 14).any()
