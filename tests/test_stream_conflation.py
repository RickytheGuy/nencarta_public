"""
Moving a source stream network onto a DEM's D8 flow paths (nencarta.tasks.stream_conflation).

Every test builds a small D8 raster whose stream paths are known cell by cell -- hillslopes drain into the
nearest stream cell -- and a source network drawn near them, and checks where each source reach lands. Whatever
else a test checks, the result must keep FLDPLN's contract: each reach is one unbroken D8 path, its downstream
reach is the one its last cell drains into, and the stream raster is exactly those cells.
"""
from collections import deque

import numpy as np
import pandas as pd
import pytest
import geopandas as gpd
from osgeo import osr
from shapely.geometry import LineString, MultiLineString

from nencarta.tasks.move_streams import _create_stream_info_table
from nencarta.tasks.stream_conflation import (
    D8_OFFSETS,
    ConflationSettings,
    FlowNetwork,
    conflate_to_flow_directions,
)

CELL = 30.0
X0, Y0 = 500_000.0, 4_500_000.0
GT = (X0, CELL, 0.0, Y0, 0.0, -CELL)
CODE = {offset: code for code, offset in D8_OFFSETS.items()}

# A threshold low enough that every cell a stream passes through is a channel cell, but most hillslope cells
# are not: a cell is 0.0009 km2, so this is about six cells' worth of drainage.
SETTINGS = ConflationSettings(min_channel_area_km2=0.005)


def utm_wkt() -> str:
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32615)
    return srs.ExportToWkt()


def centre(cell, offset=(0.0, 0.0)):
    """Map coordinates of a (row, col) cell's centre, shifted by ``offset`` (rows, cols) in cells."""
    r, c = cell
    return X0 + (c + 0.5 + offset[1]) * CELL, Y0 - (r + 0.5 + offset[0]) * CELL


def straight(r0, c0, r1, c1):
    """The cells of a straight D8 path from (r0, c0) to (r1, c1), both included."""
    cells = [(r0, c0)]
    r, c = r0, c0
    while (r, c) != (r1, c1):
        r += (r1 > r) - (r1 < r)
        c += (c1 > c) - (c1 < c)
        cells.append((r, c))
    return cells


def flowdir_for(shape, paths, exits=()):
    """
    A D8 raster where each path in ``paths`` flows from cell to cell, and every other cell drains, by the
    shortest route, into the nearest path cell. A path ending on another path's cell joins it there. The last
    cell of a path that ends nowhere has no outflow, unless it is in ``exits``, a {cell: (d_row, d_col)} of
    cells that drain off the edge of the grid.
    """
    flowdir = np.zeros(shape, np.int16)
    on_path = np.zeros(shape, bool)
    for path in paths:
        for (r0, c0), (r1, c1) in zip(path[:-1], path[1:]):
            flowdir[r0, c0] = CODE[(r1 - r0, c1 - c0)]
        for cell in path:
            on_path[cell] = True
    for cell, offset in dict(exits).items():
        flowdir[cell] = CODE[offset]
    queue = deque(zip(*np.nonzero(on_path)))
    seen = on_path.copy()
    while queue:
        r, c = queue.popleft()
        for (dr, dc) in CODE:
            nr, nc = r - dr, c - dc  # the neighbour that would drain into (r, c) with offset (dr, dc)
            if 0 <= nr < shape[0] and 0 <= nc < shape[1] and not seen[nr, nc]:
                seen[nr, nc] = True
                flowdir[nr, nc] = CODE[(dr, dc)]
                queue.append((nr, nc))
    return flowdir


def line(cells, offset=(0.0, 0.0), reverse=False):
    pts = [centre(cell, offset) for cell in cells]
    return LineString(pts[::-1] if reverse else pts)


def source(rows, crs="EPSG:32615"):
    """rows: (id, downstream id, geometry, order)."""
    return gpd.GeoDataFrame(
        {"LINKNO": [r[0] for r in rows], "DSLINKNO": [r[1] for r in rows],
         "order": [r[3] if len(r) > 3 else 1 for r in rows], "extra": [f"attr{r[0]}" for r in rows]},
        geometry=[r[2] for r in rows], crs=crs)


def conflate(src, flowdir, settings=SETTINGS, **kwargs):
    return conflate_to_flow_directions(src, flowdir, GT, utm_wkt(), "LINKNO", "DSLINKNO", order_col="order",
                                       settings=settings, **kwargs)


def cells_of(result, reach_id):
    """The raster cells carrying ``reach_id``, from upstream to downstream along the reach's line."""
    row = result.streams.set_index("LINKNO").loc[reach_id]
    out = []
    for x, y in row.geometry.coords:
        out.append((int((Y0 - y) // CELL), int((x - X0) // CELL)))
    return out


def ds_of(result, reach_id):
    return int(result.streams.set_index("LINKNO").loc[reach_id, "DSLINKNO"])


def assert_fldpln_contract(result, flowdir, tmp_path):
    """Every reach is one D8 path, drains into the reach it says it does, and the raster is exactly its cells."""
    raster = result.raster
    nrows, ncols = raster.shape
    labelled = set()
    for _, row in result.streams.iterrows():
        rid = int(row.LINKNO)
        cells = [(int((Y0 - y) // CELL), int((x - X0) // CELL)) for x, y in row.geometry.coords]
        assert len(cells) >= 2
        for a, b in zip(cells[:-1], cells[1:]):
            d_row, d_col = D8_OFFSETS[int(flowdir[a])]
            assert (a[0] + d_row, a[1] + d_col) == b, f"reach {rid} leaves its D8 path at {a}"
        for cell in cells:
            assert raster[cell] == rid
            labelled.add(cell)
        code = int(flowdir[cells[-1]])
        below = None
        if code in D8_OFFSETS:
            d_row, d_col = D8_OFFSETS[code]
            r, c = cells[-1][0] + d_row, cells[-1][1] + d_col
            if 0 <= r < nrows and 0 <= c < ncols:
                below = int(raster[r, c])
        expected_ds = below if below not in (None, 0, rid) else -1
        assert int(row.DSLINKNO) == expected_ds, f"reach {rid} says it drains into {row.DSLINKNO}, D8 says {expected_ds}"
    assert set(zip(*np.nonzero(raster))) == labelled

    # Curve2Flood sorts the network by topological_order to chain reaches: upstream reaches have to come first
    order = result.streams.set_index("LINKNO")["topological_order"]
    for rid, ds in zip(result.streams.LINKNO, result.streams.DSLINKNO):
        if ds != -1:
            assert order[rid] < order[ds]

    # And FLDPLN's walk of each reach covers exactly its cells
    flowdir_file = tmp_path / "flowdir.tif"
    if not flowdir_file.exists():
        from osgeo import gdal
        ds = gdal.GetDriverByName("GTiff").Create(str(flowdir_file), ncols, nrows, 1, gdal.GDT_Int16)
        ds.WriteArray(flowdir)
        ds = None
    table_file = tmp_path / "stream_info.parquet"
    _create_stream_info_table(table_file, raster, flowdir_file)
    table = pd.read_parquet(table_file).set_index("stream_id")
    for rid, n_cells in result.streams.set_index("LINKNO").geometry.apply(lambda g: len(g.coords)).items():
        assert table.at[rid, "length"] == n_cells


# ---------------------------------------------------------------------------------------------------------------
# One reach
# ---------------------------------------------------------------------------------------------------------------

@pytest.fixture
def river_east():
    """A river flowing east along row 10 of a 21 x 40 grid and off its east edge."""
    path = straight(10, 0, 10, 39)
    return path, flowdir_for((21, 40), [path], exits={(10, 39): (0, 1)})


def test_a_reach_on_its_d8_path_keeps_every_cell(river_east, workspace_tmp):
    path, flowdir = river_east
    result = conflate(source([(5, -1, line(path))]), flowdir)

    assert cells_of(result, 5) == path
    assert ds_of(result, 5) == -1
    assert result.streams.loc[0, "extra"] == "attr5"  # source attributes ride along
    assert_fldpln_contract(result, flowdir, workspace_tmp)


def test_a_line_drawn_downstream_to_upstream_is_turned_around(river_east, workspace_tmp):
    # TDX-Hydro and GEOGLOWS draw their lines from the outlet up
    path, flowdir = river_east
    result = conflate(source([(5, -1, line(path, reverse=True))]), flowdir)

    assert cells_of(result, 5) == path
    assert_fldpln_contract(result, flowdir, workspace_tmp)


@pytest.mark.parametrize("offset", [(1.0, 0.0), (-1.4, 0.0), (0.6, 0.3)])
def test_a_line_beside_the_channel_moves_onto_it(river_east, offset, workspace_tmp):
    path, flowdir = river_east
    result = conflate(source([(5, -1, line(path, offset=offset))]), flowdir)

    cells = cells_of(result, 5)
    assert all(r == 10 for r, _ in cells)
    assert len(cells) >= len(path) - 3
    assert_fldpln_contract(result, flowdir, workspace_tmp)


def test_a_line_far_from_any_channel_is_left_out(river_east):
    path, flowdir = river_east
    far = straight(0, 0, 0, 39)  # ten cells north, on a ridge of hillslope cells
    result = conflate(source([(5, -1, line(path)), (6, -1, line(far))]), flowdir)

    assert set(result.streams.LINKNO) == {5}
    assert not (result.raster == 6).any()


def test_a_reach_running_against_the_flow_is_left_out(river_east):
    # The D8 path runs east; a reach whose network says it runs west can't be on it
    path, flowdir = river_east
    upstream_end = straight(10, 0, 10, 19)
    result = conflate(source([(1, 2, line(path[20:])), (2, -1, line(upstream_end[::-1]))]), flowdir)

    assert 2 not in set(result.streams.LINKNO)


def test_a_geographic_grid(workspace_tmp):
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    arcsec = 1 / 3600
    gt = (-95.0, arcsec, 0.0, 41.0, 0.0, -arcsec)
    path = straight(10, 0, 10, 39)
    flowdir = flowdir_for((21, 40), [path], exits={(10, 39): (0, 1)})
    pts = [(gt[0] + (c + 0.5) * arcsec, gt[3] - (r + 0.5) * arcsec) for r, c in path]
    src = gpd.GeoDataFrame({"LINKNO": [5], "DSLINKNO": [-1]}, geometry=[LineString(pts)], crs="EPSG:4326")

    result = conflate_to_flow_directions(src, flowdir, gt, srs.ExportToWkt(), "LINKNO", "DSLINKNO",
                                         settings=ConflationSettings(min_channel_area_km2=0.003))

    assert len(result.streams) == 1
    assert int(result.raster[10, :].astype(bool).sum()) == 40


# ---------------------------------------------------------------------------------------------------------------
# Where reaches meet
# ---------------------------------------------------------------------------------------------------------------

def test_a_break_without_a_confluence_splits_the_d8_path(river_east, workspace_tmp):
    # Most of a network's reach breaks are not confluences (GEOGLOWS merges its first-order streams away)
    path, flowdir = river_east
    upper, lower = path[:21], path[20:]
    result = conflate(source([(1, 2, line(upper)), (2, -1, line(lower))]), flowdir)

    up, down = cells_of(result, 1), cells_of(result, 2)
    assert up + down == path
    assert abs(len(up) - 20) <= 1
    assert ds_of(result, 1) == 2
    assert_fldpln_contract(result, flowdir, workspace_tmp)


@pytest.fixture
def confluence():
    """
    A main stem down column 20 of a 40 x 41 grid, off the south edge, and a tributary from the north-east that
    the source network joins to it at row 15. The tributary's D8 path follows its source line, but the DEM keeps
    it beside the main stem, one column east, until row ``join``.
    """
    def build(join):
        main = straight(0, 20, 39, 20)
        tributary = straight(3, 32, 14, 21) + [(row, 21) for row in range(15, join)] + [(join, 20)]
        flowdir = flowdir_for((40, 41), [main, tributary], exits={(39, 20): (1, 0)})
        src = source([
            (1, 3, line(main[:16]), 2),     # main stem down to the source confluence
            (2, 3, line(straight(3, 32, 15, 20)), 1),
            (3, -1, line(main[15:]), 2),    # main stem below it
        ])
        return src, flowdir, main, tributary
    return build


def test_a_confluence_where_the_source_network_has_it(confluence, workspace_tmp):
    src, flowdir, main, tributary = confluence(join=15)
    result = conflate(src, flowdir)

    assert ds_of(result, 1) == 3 and ds_of(result, 2) == 3
    assert cells_of(result, 3)[0] == (15, 20)
    assert cells_of(result, 2) == tributary[:-1]
    assert_fldpln_contract(result, flowdir, workspace_tmp)


def test_a_tributary_the_dem_joins_further_down_still_joins_its_river(confluence, workspace_tmp):
    src, flowdir, main, tributary = confluence(join=22)
    result = conflate(src, flowdir)

    # It follows its own D8 path down to where that path really meets the river
    assert cells_of(result, 2)[-1] == tributary[-2]
    assert ds_of(result, 2) == 3 and ds_of(result, 1) == 3
    assert_fldpln_contract(result, flowdir, workspace_tmp)


def test_two_reaches_the_dem_joins_early_start_their_river_where_they_meet(workspace_tmp):
    # The D8 paths of the main stem and its tributary meet at row 12, three cells above the source confluence:
    # below the meeting point both flows run together, so that is where the reach they join has to start
    main = straight(0, 20, 39, 20)
    tributary = straight(3, 29, 12, 20)
    flowdir = flowdir_for((40, 41), [main, tributary], exits={(39, 20): (1, 0)})
    src = source([
        (1, 3, line(main[:16]), 2),
        (2, 3, line(straight(3, 29, 9, 23) + straight(10, 23, 15, 20)), 1),  # drawn into row 15
        (3, -1, line(main[15:]), 2),
    ])
    result = conflate(src, flowdir)

    assert cells_of(result, 3)[0] == (12, 20)
    assert cells_of(result, 1)[-1] == (11, 20)
    assert ds_of(result, 1) == 3 and ds_of(result, 2) == 3
    assert_fldpln_contract(result, flowdir, workspace_tmp)


def test_a_tributary_the_dem_routes_alongside_the_river_follows_its_own_path(workspace_tmp):
    # The DEM keeps the tributary in its own channel beside the main stem for ten cells before the two meet
    main = straight(0, 20, 39, 20)
    tributary = straight(5, 30, 12, 23) + straight(13, 23, 22, 23) + straight(23, 22, 24, 21)
    tributary = tributary + [(25, 20)]
    flowdir = flowdir_for((40, 41), [main, tributary], exits={(39, 20): (1, 0)})
    src = source([
        (1, 3, line(main[:13]), 2),
        (2, 3, line(straight(5, 30, 12, 20)), 1),  # the source joins it at row 12
        (3, -1, line(main[12:]), 2),
    ])
    result = conflate(src, flowdir)

    assert cells_of(result, 2)[-1] == (24, 21)
    assert ds_of(result, 2) == 3
    assert_fldpln_contract(result, flowdir, workspace_tmp)


def test_a_reach_that_cannot_be_matched_is_bridged_by_the_reach_above_it(river_east, workspace_tmp):
    path, flowdir = river_east
    nonsense = LineString([centre((0, 20)), centre((0, 25))])  # far from the river
    src = source([(1, 2, line(path[:15])), (2, 3, nonsense), (3, -1, line(path[25:]))])
    result = conflate(src, flowdir)

    assert set(result.streams.LINKNO) == {1, 3}
    assert cells_of(result, 1) + cells_of(result, 3) == path
    assert ds_of(result, 1) == 3
    assert_fldpln_contract(result, flowdir, workspace_tmp)


def test_two_source_lines_on_one_channel_never_share_a_cell(river_east, workspace_tmp):
    # Two unrelated source reaches drawn down the same valley: the larger keeps the channel
    path, flowdir = river_east
    src = source([(1, -1, line(path, offset=(0.4, 0)), 3), (2, -1, line(path[5:30], offset=(-0.4, 0)), 1)])
    result = conflate(src, flowdir)

    assert cells_of(result, 1) == path
    assert not (result.raster == 2).any()
    assert_fldpln_contract(result, flowdir, workspace_tmp)


# ---------------------------------------------------------------------------------------------------------------
# Endorheic basins, lakes and the edge of the domain
# ---------------------------------------------------------------------------------------------------------------

def test_an_endorheic_reach_ends_where_its_source_line_ends(workspace_tmp):
    # A filled DEM spills the closed basin at (20, 15) over its rim and on into the river along row 30; the
    # source network knows the basin is closed. Its reach must stop at the sink, and nothing may carry its id over
    # the spill path.
    river = straight(30, 0, 30, 39)
    basin = straight(5, 15, 20, 15)
    spill = straight(20, 15, 29, 15)
    flowdir = flowdir_for((40, 40), [river, basin + spill[1:] + [(30, 15)]], exits={(30, 39): (0, 1)})
    src = source([(1, -1, line(basin), 1), (7, -1, line(river), 2)])
    result = conflate(src, flowdir)

    basin_cells = cells_of(result, 1)
    assert basin_cells[-1][0] in (19, 20, 21)
    assert ds_of(result, 1) == -1
    assert not result.raster[23:30, 15].any()
    assert cells_of(result, 7) == river
    assert_fldpln_contract(result, flowdir, workspace_tmp)


@pytest.fixture
def river_through_lake():
    """A river flowing east along row 10 through a lake covering columns 15 to 24."""
    path = straight(10, 0, 10, 39)
    flowdir = flowdir_for((21, 40), [path], exits={(10, 39): (0, 1)})
    lake = np.zeros((21, 40), bool)
    lake[6:15, 15:25] = True
    src = source([
        (1, 2, line(path[:18])),     # flows into the lake
        (2, 3, line(path[17:23])),   # inside it
        (3, -1, line(path[22:])),    # leaves it
    ])
    return path, flowdir, lake, src


def test_a_lake_cuts_the_network_at_its_shores(river_through_lake, workspace_tmp):
    path, flowdir, lake, src = river_through_lake
    result = conflate(src, flowdir, lake_mask=lake)

    assert set(result.streams.LINKNO) == {1, 3}
    assert cells_of(result, 1)[-1] == (10, 14)
    assert cells_of(result, 3)[0] == (10, 25)
    assert not (result.raster.astype(bool) & lake).any()
    assert ds_of(result, 1) == -1
    assert_fldpln_contract(result, flowdir, workspace_tmp)


def test_without_a_lake_mask_the_lake_reach_is_kept(river_through_lake, workspace_tmp):
    path, flowdir, lake, src = river_through_lake
    result = conflate(src, flowdir)

    assert set(result.streams.LINKNO) == {1, 2, 3}
    assert cells_of(result, 1) + cells_of(result, 2) + cells_of(result, 3) == path


def test_a_river_flowing_in_across_the_edge_is_matched_from_the_edge(workspace_tmp):
    # With a channel threshold larger than the strip of domain the river has crossed, the river's first cells
    # drain too little to count as channel -- unless its inflow is counted
    path = straight(10, 0, 10, 39)
    flowdir = flowdir_for((21, 40), [path], exits={(10, 39): (0, 1)})
    entering = LineString([centre((10, -8)), centre((10, 39))])  # starts outside the grid
    settings = ConflationSettings(min_channel_area_km2=0.2)
    result = conflate(source([(5, -1, entering)]), flowdir, settings=settings)

    assert cells_of(result, 5)[0] == (10, 0)
    assert_fldpln_contract(result, flowdir, workspace_tmp)


def test_nodata_stops_the_network(river_east, workspace_tmp):
    path, flowdir = river_east
    valid = np.ones(flowdir.shape, bool)
    valid[:, 30:] = False  # ocean
    result = conflate(source([(5, -1, line(path))]), flowdir, valid_mask=valid)

    assert cells_of(result, 5)[-1] == (10, 29)
    assert not result.raster[:, 30:].any()


# ---------------------------------------------------------------------------------------------------------------
# Messy inputs
# ---------------------------------------------------------------------------------------------------------------

def test_messy_source_networks_are_tolerated(river_east, workspace_tmp):
    path, flowdir = river_east
    rows = [
        (1, 2, line(path[:10])),
        (2, 1, line(path[9:20])),                          # a cycle with reach 1
        (3, 3, line(path[19:30])),                         # drains into itself
        (4, 999, line(path[29:])),                         # drains into a reach that isn't here
        (4, -1, line(path[29:])),                          # listed twice
        (0, -1, line(path[:5])),                           # an id that can't label a raster cell
        (-7, -1, line(path[:5])),
        (8, -1, None),                                     # no geometry
        (9, -1, LineString()),                             # empty geometry
        (10, -1, MultiLineString([line(path[:3]), line(path[30:33])])),  # in two pieces
    ]
    src = gpd.GeoDataFrame({"LINKNO": [r[0] for r in rows], "DSLINKNO": [r[1] for r in rows], "order": 1},
                           geometry=[r[2] for r in rows], crs="EPSG:32615")
    result = conflate(src, flowdir)

    assert {1, 2, 3, 4} <= set(result.streams.LINKNO)
    assert result.streams.LINKNO.is_unique
    assert not (result.raster < 0).any()
    assert_fldpln_contract(result, flowdir, workspace_tmp)


@pytest.mark.parametrize("has_downstream", [True, False])
def test_a_y_shaped_reach_is_matched_along_its_trunk_and_longer_branch(has_downstream, workspace_tmp):
    # TDX-Hydro folds first-order streams into the reach below, leaving reaches drawn as three pieces that
    # line_merge can't join: two headwater branches and the trunk below where they meet. The trunk is what has
    # to be matched -- it carries the flow -- and it is not the longest piece.
    trunk = straight(20, 20, 20, 39)
    long_branch = straight(2, 2, 20, 20)
    short_branch = straight(30, 10, 20, 20)
    main = long_branch + trunk[1:]
    paths = [main, short_branch]
    exits = {(20, 39): (0, 1)}
    rows = []
    y_reach = MultiLineString([line(trunk, reverse=True), line(long_branch, reverse=True),
                               line(short_branch, reverse=True)])
    if has_downstream:
        below = straight(20, 39, 20, 49)
        paths = [long_branch + trunk[1:] + below[1:], short_branch]
        exits = {(20, 49): (0, 1)}
        rows.append((2, -1, line(below), 1))
    rows.insert(0, (1, 2 if has_downstream else -1, y_reach, 1))
    flowdir = flowdir_for((41, 50 if has_downstream else 40), paths, exits=exits)
    result = conflate(source(rows), flowdir)

    cells = cells_of(result, 1)
    assert cells[0] == long_branch[0]
    if has_downstream:
        assert cells[-1] == trunk[-2]  # the reach below starts at the trunk's last cell
        assert ds_of(result, 1) == 2
    else:
        assert cells[-1] == trunk[-1]
    assert_fldpln_contract(result, flowdir, workspace_tmp)


def test_a_flow_cycle_in_the_pointer_raster_does_not_hang(river_east, workspace_tmp):
    path, flowdir = river_east
    flowdir = flowdir.copy()
    flowdir[0, 0] = CODE[(0, 1)]
    flowdir[0, 1] = CODE[(0, -1)]
    result = conflate(source([(5, -1, line(path))]), flowdir)

    assert cells_of(result, 5) == path


def test_ids_too_large_for_the_stream_raster_raise(river_east):
    path, flowdir = river_east
    with pytest.raises(ValueError, match="32-bit"):
        conflate(source([(2**31, -1, line(path))]), flowdir)


def test_an_empty_source_network(river_east):
    path, flowdir = river_east
    empty = gpd.GeoDataFrame({"LINKNO": [], "DSLINKNO": [], "order": []}, geometry=[], crs="EPSG:32615")
    result = conflate(empty, flowdir)

    assert result.streams.empty
    assert not result.raster.any()


def test_a_reach_too_short_to_match_hands_its_cells_to_the_reach_above(river_east, workspace_tmp):
    path, flowdir = river_east
    src = source([(1, 2, line(path[:20])), (2, 3, line(path[19:21])), (3, -1, line(path[20:]))])
    result = conflate(src, flowdir)

    assert 2 not in set(result.streams.LINKNO)
    assert cells_of(result, 1) + cells_of(result, 3) == path
    assert ds_of(result, 1) == 3
    assert_fldpln_contract(result, flowdir, workspace_tmp)


def test_a_large_tree_keeps_the_contract(workspace_tmp):
    # A trunk with eight side tributaries, each source line nudged off its channel
    rng = np.random.default_rng(3)
    trunk = straight(30, 0, 30, 79)
    paths, rows = [trunk], []
    breaks = list(range(0, 80, 10)) + [80]
    for k, (a, b) in enumerate(zip(breaks[:-1], breaks[1:]), start=1):
        rows.append((100 + k, 100 + k + 1 if b < 80 else -1, line(trunk[a:b + 1] if b < 80 else trunk[a:]), 2))
    for k, col in enumerate(range(5, 80, 10), start=1):
        side = 1 if k % 2 else -1
        trib = straight(30 + side * 12, col, 30 + side, col) + [(30, col)]
        paths.append(trib)
        reach = 100 + 1 + col // 10
        rows.append((k, reach, line(trib[:-1], offset=tuple(rng.uniform(-0.6, 0.6, 2))), 1))
    flowdir = flowdir_for((61, 80), paths, exits={(30, 79): (0, 1)})
    result = conflate(source(rows), flowdir)

    assert len(result.streams) == len(rows)
    assert_fldpln_contract(result, flowdir, workspace_tmp)


def test_a_flow_network_can_be_reused(river_east):
    path, flowdir = river_east
    flow = FlowNetwork(flowdir, GT, utm_wkt())
    first = conflate(source([(5, -1, line(path))]), flowdir, flow=flow)
    second = conflate(source([(5, -1, line(path))]), flowdir, flow=flow)

    assert np.array_equal(first.raster, second.raster)
