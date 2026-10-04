"""
Move a source stream network onto a DEM's D8 flow paths.

FLDPLN, ARC and Curve2Flood all read the stream network off the flow direction raster: FLDPLN walks each reach
down the D8 pointers, and a reach that leaves its D8 path mid-way stops being one reach. A source network
(TDX-Hydro, GEOGLOWS, NHD...) was drawn from some other DEM, so its lines drift on and off this DEM's flow
paths. This module moves every source reach onto the D8 path it follows, keeping its id so the flows keyed by
that id still apply.

It works from the source network down, not from a DEM-derived network up:

1. Every source reach is sampled about once a cell, from its upstream end down, and matched onto the D8 network
   as a hidden Markov model (map matching, Newson & Krumm 2009). The hidden states at a sample are the channel
   cells near it; a state can only be followed by one downstream of it on the D8 network; and a transition costs
   however much farther, or shorter, the flow path runs between the two cells than the source line runs between
   the two samples. The best local alignment is kept, so a reach that only partly follows the DEM keeps the part
   that does, and a reach that doesn't follow it at all is left out rather than forced somewhere wrong.
2. The matched reaches are painted onto the D8 network in flow order. A reach runs from its matched start down
   to where the next reach starts. Where the flow paths of two reaches meet somewhere the source network has
   no confluence, the source network decides who carries on: a reach downstream of the others in the source
   network wins, two reaches that the source network joins into a third start that third reach where they
   meet, and otherwise the larger river wins. A reach that meets no other reach ends where its match ended.

Because every reach is matched on its own, a reach the DEM disagrees with costs only that reach, never the
basin above it. The consequences, which the alternative of labelling a DEM-derived network did not have:

- Reaches split exactly where the source network splits them, including the many breaks that are not
  confluences, without having to cut a DEM-derived network first.
- An endorheic source reach ends where its source line ends, even though a filled DEM routes its flow onwards
  over the basin's spill point: no reach is downstream of it, so nothing carries its id past its own match.
- A river leaving the domain, or flowing out through the wrong side of a flat reservoir, is still matched up to
  where it leaves.
- Lakes can be cut out: cells inside a lake are never matched and stop a reach that flows into them.

The stream raster is the painted D8 network itself, so every reach is one unbroken D8 path.
"""
from __future__ import annotations

import heapq
import math
from dataclasses import dataclass

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
from numba import njit
from osgeo import osr

from nencarta.logger import LOG

# WhiteboxTools D8 pointer codes, mapped to the (row, col) offset of the cell each
# code drains into. Curve2Flood's FLDPLN spreader walks the stream info table with
# this same encoding, so the table has to be built against it.
D8_OFFSETS = {
    64: (-1, -1), 128: (-1, 0), 1: (-1, 1),
    32: (0, -1),                2: (0, 1),
    16: (1, -1),    8: (1, 0),  4: (1, 1),
}


def d8_offset_tables() -> tuple[np.ndarray, np.ndarray]:
    """Row/column offset per D8 code. An offset of (0, 0) means 'no outflow'."""
    row_offsets = np.zeros(256, np.int64)
    col_offsets = np.zeros(256, np.int64)
    for code, (d_row, d_col) in D8_OFFSETS.items():
        row_offsets[code] = d_row
        col_offsets[code] = d_col
    return row_offsets, col_offsets


@dataclass
class ConflationSettings:
    """
    How the source network is matched onto the D8 network. Distances are in DEM cells, except those named in metres.

    The defaults suit a ~30 m DEM with the source streams burned into it, where the D8 paths and the source lines
    rarely disagree by more than a cell or two.
    """
    # A cell is a channel cell, and so a candidate for a source line to match onto, when this much area drains
    # through it. Well below where any source network starts its streams, so that small source headwaters still
    # find a channel to match.
    min_channel_area_km2: float = 0.25
    # How far from a source sample a channel cell may be and still be a candidate, and how many of the nearest
    # candidates are kept. The radius grows with the river (see search_radius_per_sqrt_km2).
    search_radius_cells: float = 5.0
    max_candidates: int = 16
    # Wide rivers' centrelines and thalwegs can be far apart: the radius grows by this many metres per square
    # root of drainage area (km2), up to max_search_radius_m. 1.0 gives ~300 m at 100 000 km2.
    search_radius_per_sqrt_km2: float = 1.0
    max_search_radius_m: float = 600.0
    # Sample spacing along a source line, in cells.
    sample_spacing_cells: float = 1.0
    # A sample's fit to a cell is 1 - (distance / sigma)^2 / 2: positive within ~1.4 sigma.
    sigma_cells: float = 1.5
    # A transition costs (|D8 path length - source line length| - slack) / beta between consecutive samples.
    beta_cells: float = 1.0
    slack_cells: float = 0.5
    # How many samples a match may skip in a row (a short stretch where the source line strays off the channel),
    # and what each skipped sample costs on top of losing its fit.
    max_gap_samples: int = 6
    gap_penalty: float = 0.25
    # A match must cover at least this many samples, and this fraction of the samples that could have been matched
    # (inside the domain and outside lakes), to be kept. On real tiles the matches below a quarter are a few
    # hundred metres of a reach that hardly follows the DEM, nearly all in the buffer along the domain's edge.
    min_matched_samples: int = 3
    min_matched_fraction: float = 0.25
    # A reach whose match ends short of the reach the source network sends it to may carry on down its D8 path to
    # reach it, for at most this far. Covers the gap left by a short reach that could not be matched, and a
    # tributary the DEM routes alongside the main stem before the two meet.
    max_extension_m: float = 3000.0
    # A reach whose match ends this close to some other matched reach is joined to it, whatever the source
    # network says about the two.
    snap_cells: float = 3.0
    # Reaches covering fewer cells than this are folded into the reach upstream of them.
    min_reach_cells: int = 2


# ---------------------------------------------------------------------------------------------------------------
# The D8 network
# ---------------------------------------------------------------------------------------------------------------

def cell_sizes_m(geotransform: tuple, projection: str, nrows: int) -> tuple[np.ndarray, float]:
    """
    Width of a cell in metres on each row, and the height of a cell in metres.

    For a geographic CRS the width shrinks with the cosine of the latitude, so it is given per row; the height
    barely changes over a tile and is taken at its centre. Rotated grids are not supported.
    """
    if geotransform[2] != 0 or geotransform[4] != 0:
        raise ValueError(f"Rotated grids are not supported: geotransform {geotransform}.")
    srs = osr.SpatialReference()
    srs.ImportFromWkt(projection)
    if srs.IsGeographic():
        lat = np.radians(geotransform[3] + (np.arange(nrows) + 0.5) * geotransform[5])
        metres_per_degree_lat = 111132.92 - 559.82 * np.cos(2 * lat) + 1.175 * np.cos(4 * lat)
        metres_per_degree_lon = 111412.84 * np.cos(lat) - 93.5 * np.cos(3 * lat)
        dx = np.abs(geotransform[1]) * metres_per_degree_lon
        dy = float(np.abs(geotransform[5]) * metres_per_degree_lat[nrows // 2])
        return dx.astype(np.float64), dy
    metres_per_unit = srs.GetLinearUnits() if srs.IsProjected() else 1.0
    dx = np.full(nrows, abs(geotransform[1]) * metres_per_unit, np.float64)
    return dx, float(abs(geotransform[5]) * metres_per_unit)


@njit(cache=True)
def _downstream_cells(flowdir: np.ndarray, valid: np.ndarray, row_offsets: np.ndarray,
                      col_offsets: np.ndarray) -> np.ndarray:
    """The flat index of the cell each cell drains into, or -1 where its flow stops (a pit, nodata, the edge)."""
    nrows, ncols = flowdir.shape
    down = np.full(nrows * ncols, -1, np.int64)
    for row in range(nrows):
        for col in range(ncols):
            if not valid[row, col]:
                continue
            code = flowdir[row, col]
            if code < 0 or code > 255:
                continue
            d_row = row_offsets[code]
            d_col = col_offsets[code]
            if d_row == 0 and d_col == 0:
                continue
            to_row = row + d_row
            to_col = col + d_col
            if to_row < 0 or to_row >= nrows or to_col < 0 or to_col >= ncols or not valid[to_row, to_col]:
                continue
            down[row * ncols + col] = to_row * ncols + to_col
    return down


@njit(cache=True)
def _accumulate(down: np.ndarray, cell_area_by_row: np.ndarray, ncols: int) -> tuple[np.ndarray, int]:
    """
    Drainage area of every cell, in flow order. Cells on a flow cycle, which a pointer raster made from a filled
    DEM never has, are cut loose from their downstream cell (in place) so that everything else stays a forest.
    Returns the areas and the number of cells cut.
    """
    n = down.size
    pending = np.zeros(n, np.int32)
    for i in range(n):
        j = down[i]
        if j >= 0:
            pending[j] += 1
    area = np.empty(n, np.float32)
    for i in range(n):
        area[i] = cell_area_by_row[i // ncols]
    queue = np.empty(n, np.int64)
    head = 0
    tail = 0
    for i in range(n):
        if pending[i] == 0:
            queue[tail] = i
            tail += 1
    while head < tail:
        i = queue[head]
        head += 1
        j = down[i]
        if j >= 0:
            area[j] += area[i]
            pending[j] -= 1
            if pending[j] == 0:
                queue[tail] = j
                tail += 1
    cut = 0
    if tail < n:
        for i in range(n):
            if pending[i] > 0 and down[i] >= 0:
                down[i] = -1
                cut += 1
    return area, cut


@njit(cache=True)
def _euler_tour(parent: np.ndarray, step: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Number a forest given by each node's parent (-1 at a root) so that ``b`` is an ancestor of ``a``, or ``a``
    itself, exactly when ``tin[b] <= tin[a] <= tout[b]``, and give each node its summed ``step`` to its root.

    With parents pointing downstream that reads "``b`` is downstream of ``a``". Iterative: numba's cache cannot
    reload a recursive function. Nodes that never reach a root keep ``tin == -1``.
    """
    m = parent.size
    first_child = np.zeros(m + 1, np.int64)
    for k in range(m):
        p = parent[k]
        if p >= 0:
            first_child[p + 1] += 1
    for k in range(m):
        first_child[k + 1] += first_child[k]
    children = np.empty(max(first_child[m], 1), np.int64)
    fill = first_child[:m].copy()
    for k in range(m):
        p = parent[k]
        if p >= 0:
            children[fill[p]] = k
            fill[p] += 1

    tin = np.full(m, -1, np.int64)
    tout = np.full(m, -1, np.int64)
    to_root = np.zeros(m, np.float64)
    stack = np.empty(m, np.int64)
    next_child = np.empty(m, np.int64)
    clock = 0
    for root in range(m):
        if parent[root] >= 0:
            continue
        stack[0] = root
        next_child[root] = first_child[root]
        tin[root] = clock
        clock += 1
        top = 0
        while top >= 0:
            k = stack[top]
            if next_child[k] < first_child[k + 1]:
                child = children[next_child[k]]
                next_child[k] += 1
                to_root[child] = to_root[k] + step[child]
                tin[child] = clock
                clock += 1
                top += 1
                stack[top] = child
                next_child[child] = first_child[child]
            else:
                tout[k] = clock - 1
                top -= 1
    return tin, tout, to_root


@njit(cache=True)
def _step_lengths(cells: np.ndarray, down: np.ndarray, ncols: int, dx_m: np.ndarray, dy_m: float) -> np.ndarray:
    """Metres from each of ``cells`` to the cell it drains into (0 where its flow stops)."""
    out = np.zeros(cells.size, np.float64)
    for k in range(cells.size):
        i = cells[k]
        j = down[i]
        if j < 0:
            continue
        d_row = j // ncols - i // ncols
        d_col = j % ncols - i % ncols
        out[k] = math.sqrt((d_col * dx_m[i // ncols]) ** 2 + (d_row * dy_m) ** 2)
    return out


@njit(cache=True)
def _add_inflow(down: np.ndarray, area: np.ndarray, cells: np.ndarray, amount: float) -> None:
    """Add ``amount`` to the drainage area of each of ``cells`` and of every cell downstream of it, in place."""
    for k in range(cells.size):
        cell = cells[k]
        steps = 0
        while cell >= 0 and steps < down.size:
            area[cell] += amount
            cell = down[cell]
            steps += 1


class FlowNetwork:
    """
    A D8 flow direction raster as a forest of cells, each pointing at the cell it drains into.

    After ``build_channels``, cells draining at least a threshold area are channel cells. They are closed
    downstream -- a channel cell drains into another channel cell -- so they form a forest of their own, which is
    given an Euler tour: whether one channel cell is downstream of another, and how far along the flow path, are
    then O(1).
    """

    def __init__(self, flowdir: np.ndarray, geotransform: tuple, projection: str, valid: np.ndarray | None = None):
        flowdir = np.ascontiguousarray(flowdir)
        if flowdir.dtype.kind not in "iu":
            flowdir = flowdir.astype(np.int64)
        self.nrows, self.ncols = flowdir.shape
        self.geotransform = tuple(geotransform)
        self.projection = projection
        if valid is None:
            valid = (flowdir >= 0) & (flowdir <= 255)
        self.valid = np.ascontiguousarray(valid, dtype=np.bool_)
        self.dx_m, self.dy_m = cell_sizes_m(self.geotransform, projection, self.nrows)
        self.cell_size_m = float(math.sqrt(float(np.median(self.dx_m)) * self.dy_m))

        row_offsets, col_offsets = d8_offset_tables()
        self.down = _downstream_cells(flowdir, self.valid, row_offsets, col_offsets)
        self.area_km2, cut = _accumulate(self.down, self.dx_m * self.dy_m / 1e6, self.ncols)
        if cut:
            LOG.warning(f"The flow direction raster has {cut} cells on flow cycles; they were cut from the network.")
        self.channel_cells = np.empty(0, np.int64)

    def build_channels(self, min_channel_area_km2: float, inflow_cells: np.ndarray | None = None) -> None:
        """
        Make the cells draining at least ``min_channel_area_km2`` the channel network.

        A river flowing in across the edge of the domain drains almost nothing within it where it enters, so the
        flow path from each of ``inflow_cells`` is counted as draining the threshold area on top of its own.
        """
        area = self.area_km2
        if inflow_cells is not None and len(inflow_cells):
            area = area.copy()
            _add_inflow(self.down, area, np.asarray(inflow_cells, np.int64), np.float32(min_channel_area_km2))
        self.channel_cells = np.flatnonzero((area >= min_channel_area_km2) & self.valid.ravel())
        self.channel_id = np.full(self.down.size, -1, np.int64)
        self.channel_id[self.channel_cells] = np.arange(self.channel_cells.size)
        downstream = self.down[self.channel_cells]
        parent = np.where(downstream >= 0, self.channel_id[np.maximum(downstream, 0)], -1)
        step = _step_lengths(self.channel_cells, self.down, self.ncols, self.dx_m, self.dy_m)
        self.tin, self.tout, self.to_outlet_m = _euler_tour(parent, step)

    def cell_centres(self, cells: np.ndarray) -> np.ndarray:
        """Map coordinates of the centres of flat ``cells``, as an (n, 2) array."""
        gt = self.geotransform
        rows = cells // self.ncols + 0.5
        cols = cells % self.ncols + 0.5
        return np.column_stack((gt[0] + cols * gt[1], gt[3] + rows * gt[5]))


# ---------------------------------------------------------------------------------------------------------------
# The source network
# ---------------------------------------------------------------------------------------------------------------

@dataclass
class SourceNetwork:
    """The source reaches, indexed 0..n-1, with their downstream reach and a numbering for ancestry tests."""
    ids: np.ndarray            # int64 source id of each reach
    ds: np.ndarray             # int64 index of the reach downstream of each reach, -1 for none
    geometry: np.ndarray       # shapely LineStrings, drawn upstream to downstream
    order: np.ndarray          # float64 stream order (0 where unknown)
    area: np.ndarray           # float64 drainage area in km2 (0 where unknown)
    tin: np.ndarray            # Euler tour of the network with parents downstream: a is b or downstream of it
    tout: np.ndarray           # exactly when tin[a] <= tin[b] <= tout[a]


def _line_parts(geom) -> list:
    """The non-empty LineStrings of a (merged) line geometry."""
    if geom is None or geom.is_empty:
        return []
    if geom.geom_type == "LineString":
        return [geom]
    return [g for g in getattr(geom, "geoms", []) if g.geom_type == "LineString" and not g.is_empty and g.length > 0]


def _main_stem(parts: list, outlet_score) -> shapely.LineString:
    """
    One line through a reach drawn in several pieces: from its outlet up the longest way through the pieces,
    drawn from upstream down.

    A reach whose pieces line_merge can't join is usually Y-shaped. TDX-Hydro folds first-order streams into the
    reach below them, which leaves a quarter of the reaches at the FIM benchmark sites as two headwater branches
    and the trunk below where they meet. The trunk carries the reach's flow, so it is the part that has to be
    matched, and the longer branch carries it on upstream; the other branch can't be on the same D8 path. Pieces
    are joined where their ends touch; ``outlet_score`` ranks the ends, and the highest is the outlet.
    """
    tolerance = 1e-6 * sum(p.length for p in parts)
    ends = [(np.asarray(p.coords[0][:2]), np.asarray(p.coords[-1][:2])) for p in parts]
    nodes = []
    part_nodes = []
    for a, b in ends:
        found = []
        for xy in (a, b):
            for k, node in enumerate(nodes):
                if np.hypot(*(node - xy)) <= tolerance:
                    found.append(k)
                    break
            else:
                nodes.append(xy)
                found.append(len(nodes) - 1)
        part_nodes.append(tuple(found))
    scores = [outlet_score(shapely.Point(xy)) for xy in nodes]
    outlet = int(np.argmax(scores))

    # The farthest node from the outlet through the pieces, and the pieces on the way there
    adjacent = [[] for _ in nodes]
    for i, (a, b) in enumerate(part_nodes):
        adjacent[a].append((i, b))
        adjacent[b].append((i, a))
    reach = {outlet: (0.0, [])}
    stack = [outlet]
    while stack:
        node = stack.pop()
        length, path = reach[node]
        for i, other in adjacent[node]:
            if other not in reach:
                reach[other] = (length + parts[i].length, path + [(i, node)])
                stack.append(other)
    far = max(reach, key=lambda k: reach[k][0])
    path = reach[far][1]
    if not path:  # the outlet's piece touches no other: take the longest piece on its own
        part = max(range(len(parts)), key=lambda i: (outlet in part_nodes[i], parts[i].length))
        coords = np.asarray(parts[part].coords)[:, :2]
        return shapely.LineString(coords if part_nodes[part][1] == outlet else coords[::-1])

    # Walk back from the far end to the outlet, each piece turned to run that way
    coords = []
    for i, downstream_node in reversed(path):
        piece = np.asarray(parts[i].coords)[:, :2]
        if part_nodes[i][1] != downstream_node:
            piece = piece[::-1]
        coords.extend(piece if not coords else piece[1:])
    return shapely.LineString(coords)


def _break_cycles(ds: np.ndarray) -> int:
    """Point the reach that closes each cycle of ``ds`` at nothing. Returns how many cycles were broken."""
    n = ds.size
    state = np.zeros(n, np.int8)  # 0 unseen, 1 on the current walk, 2 done
    broken = 0
    for start in range(n):
        walk = []
        node = start
        while node >= 0 and state[node] == 0:
            state[node] = 1
            walk.append(node)
            node = ds[node]
        if node >= 0 and state[node] == 1:
            ds[walk[-1]] = -1
            broken += 1
        for visited in walk:
            state[visited] = 2
    return broken


def prepare_source_network(source_gdf: gpd.GeoDataFrame, id_col: str, ds_col: str, order_col: str | None,
                           area_col: str | None, flow: FlowNetwork) -> SourceNetwork:
    """
    Clean the source reaches into a forest whose lines run downstream.

    Reaches without a usable line or with an id that can't label a raster cell (missing, or not above 0) are
    dropped; a reach listed twice keeps its first row; a downstream id that isn't a reach here (outside the
    domain, or 0/-1 for "none") becomes "no downstream"; a reach pointing at itself or closing a cycle is cut.

    Source networks draw their lines either way (TDX-Hydro and GEOGLOWS draw them downstream to upstream), so
    each line is turned to run downstream. A line's downstream end is the one nearer the reach it drains into;
    without one, the end farther from the reaches draining into it; and for a reach with neither, the end with
    more area draining through the DEM there.
    """
    gdf = source_gdf[[c for c in dict.fromkeys([id_col, ds_col, order_col, area_col, source_gdf.geometry.name]) if c]]
    ids = pd.to_numeric(gdf[id_col], errors="coerce")
    geoms = gdf.geometry.to_numpy()
    usable = (ids.notna() & (ids > 0)).to_numpy() & ~shapely.is_missing(geoms) & ~shapely.is_empty(geoms)
    gdf = gdf[usable]
    gdf = gdf[~pd.to_numeric(gdf[id_col]).duplicated().to_numpy()]

    merged = shapely.line_merge(shapely.force_2d(gdf.geometry.to_numpy()))
    parts = [_line_parts(g) for g in merged]
    keep = np.array([len(p) > 0 for p in parts], bool)
    gdf, merged = gdf[keep], merged[keep]
    parts = [p for p, k in zip(parts, keep) if k]

    ids = pd.to_numeric(gdf[id_col]).to_numpy(np.int64)
    position = {rid: i for i, rid in enumerate(ids)}
    ds_ids = pd.to_numeric(gdf[ds_col], errors="coerce").to_numpy()
    ds = np.array([position.get(int(d), -1) if np.isfinite(d) else -1 for d in ds_ids], np.int64)
    ds[ds == np.arange(ds.size)] = -1
    broken = _break_cycles(ds)
    if broken:
        LOG.warning(f"The source network has {broken} cycles; each was cut at one reach.")

    # A reach whose line can't be merged into one is matched along its main stem (see _main_stem). Its outlet is
    # the end nearest the reach it drains into; without one, the end farthest from the reaches draining into it;
    # and with neither, the end the most drains through on the DEM.
    upstream = [[] for _ in range(ids.size)]
    for i in np.flatnonzero(ds >= 0):
        upstream[ds[i]].append(i)
    lines = np.empty(ids.size, dtype=object)
    for i, pieces in enumerate(parts):
        if len(pieces) == 1:
            lines[i] = pieces[0]
            continue
        if ds[i] >= 0:
            target = merged[ds[i]]
            lines[i] = _main_stem(pieces, lambda point: -point.distance(target))
        elif upstream[i]:
            sources = shapely.union_all(merged[upstream[i]])
            lines[i] = _main_stem(pieces, lambda point: point.distance(sources))
        else:
            lines[i] = _main_stem(pieces, lambda point: _area_near(flow, point))

    order = pd.to_numeric(gdf[order_col], errors="coerce").fillna(0).to_numpy(np.float64) \
        if order_col and order_col in gdf.columns else np.zeros(ids.size)
    area = pd.to_numeric(gdf[area_col], errors="coerce").fillna(0).to_numpy(np.float64) \
        if area_col and area_col in gdf.columns else np.zeros(ids.size)

    lines = _orient_downstream(lines, ds, flow)
    tin, tout, _ = _euler_tour(ds, np.zeros(ds.size))
    return SourceNetwork(ids, ds, lines, order, area, tin, tout)


def _orient_downstream(lines: np.ndarray, ds: np.ndarray, flow: FlowNetwork) -> np.ndarray:
    n = lines.size
    first = shapely.points([g.coords[0] for g in lines]) if n else np.array([], dtype=object)
    last = shapely.points([g.coords[-1] for g in lines]) if n else np.array([], dtype=object)
    flip = np.zeros(n, bool)
    decided = np.zeros(n, bool)

    has_ds = ds >= 0
    if has_ds.any():
        target = lines[ds[has_ds]]
        d_first = shapely.distance(first[has_ds], target)
        d_last = shapely.distance(last[has_ds], target)
        idx = np.flatnonzero(has_ds)
        flip[idx] = d_first < d_last
        decided[idx] = d_first != d_last

    upstream = [[] for _ in range(n)]
    for i in np.flatnonzero(has_ds):
        upstream[ds[i]].append(i)
    for i in np.flatnonzero(~decided):
        if not upstream[i]:
            continue
        ups = shapely.union_all(lines[upstream[i]])
        d_first = shapely.distance(first[i], ups)
        d_last = shapely.distance(last[i], ups)
        if d_first != d_last:
            flip[i] = d_last < d_first
            decided[i] = True

    # A reach with neither: the end where more drains through the DEM. Taken at the line's last vertices on the
    # DEM, since an end can lie off the grid or out at sea.
    for i in np.flatnonzero(~decided):
        on_dem = [p for p in shapely.points(np.asarray(lines[i].coords)[:, :2]) if _on_valid_cell(flow, p)]
        if len(on_dem) >= 2:
            flip[i] = _area_near(flow, on_dem[0]) > _area_near(flow, on_dem[-1])

    out = lines.copy()
    for i in np.flatnonzero(flip):
        out[i] = shapely.reverse(lines[i])
    return out


def _on_valid_cell(flow: FlowNetwork, point) -> bool:
    gt = flow.geotransform
    col = int(math.floor((point.x - gt[0]) / gt[1]))
    row = int(math.floor((point.y - gt[3]) / gt[5]))
    return 0 <= row < flow.nrows and 0 <= col < flow.ncols and bool(flow.valid[row, col])


def _area_near(flow: FlowNetwork, point, radius_cells: int = 2) -> float:
    """The largest drainage area on the DEM within a couple of cells of ``point``, or 0 if there is none."""
    gt = flow.geotransform
    col = int(math.floor((point.x - gt[0]) / gt[1]))
    row = int(math.floor((point.y - gt[3]) / gt[5]))
    r0, r1 = max(row - radius_cells, 0), min(row + radius_cells + 1, flow.nrows)
    c0, c1 = max(col - radius_cells, 0), min(col + radius_cells + 1, flow.ncols)
    if r0 >= r1 or c0 >= c1:
        return 0.0
    window = flow.area_km2.reshape(flow.nrows, flow.ncols)[r0:r1, c0:c1]
    valid = flow.valid[r0:r1, c0:c1]
    return float(window[valid].max()) if valid.any() else 0.0


# ---------------------------------------------------------------------------------------------------------------
# Matching each reach
# ---------------------------------------------------------------------------------------------------------------

@njit(cache=True)
def _resample(xs: np.ndarray, ys: np.ndarray, offsets: np.ndarray, spacing: float
              ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Points every ``spacing`` along each polyline (vertices ``offsets[i]:offsets[i+1]``, in metre-like
    coordinates), always including both ends. Returns the points, the distance of each along its line, and the
    offsets of each line's points.
    """
    n_lines = offsets.size - 1
    total = 0
    for i in range(n_lines):
        length = 0.0
        for v in range(offsets[i] + 1, offsets[i + 1]):
            length += math.hypot(xs[v] - xs[v - 1], ys[v] - ys[v - 1])
        total += int(math.floor(length / spacing)) + 2
    out_x = np.empty(total, np.float64)
    out_y = np.empty(total, np.float64)
    along = np.empty(total, np.float64)
    out_offsets = np.zeros(n_lines + 1, np.int64)
    k = 0
    for i in range(n_lines):
        start, stop = offsets[i], offsets[i + 1]
        out_x[k] = xs[start]
        out_y[k] = ys[start]
        along[k] = 0.0
        k += 1
        travelled = 0.0
        next_at = spacing
        for v in range(start + 1, stop):
            seg = math.hypot(xs[v] - xs[v - 1], ys[v] - ys[v - 1])
            while seg > 0 and next_at <= travelled + seg - 1e-9:
                t = (next_at - travelled) / seg
                out_x[k] = xs[v - 1] + t * (xs[v] - xs[v - 1])
                out_y[k] = ys[v - 1] + t * (ys[v] - ys[v - 1])
                along[k] = next_at
                k += 1
                next_at += spacing
            travelled += seg
        if travelled > along[k - 1] + 1e-9 or k - out_offsets[i] == 1:
            out_x[k] = xs[stop - 1]
            out_y[k] = ys[stop - 1]
            along[k] = travelled
            k += 1
        out_offsets[i + 1] = k
    return out_x[:k], out_y[:k], along[:k], out_offsets


@njit(cache=True)
def _find_candidates(rows_f: np.ndarray, cols_f: np.ndarray, usable: np.ndarray, radius_m: np.ndarray,
                     k_max: int, channel_id: np.ndarray, excluded: np.ndarray, nrows: int, ncols: int,
                     dx_m: np.ndarray, dy_m: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The ``k_max`` nearest channel cells (as channel ids) within ``radius_m`` of each usable sample."""
    n = rows_f.size
    cand = np.full((n, k_max), -1, np.int64)
    cand_dist = np.full((n, k_max), np.inf, np.float64)
    count = np.zeros(n, np.int64)
    for s in range(n):
        if not usable[s]:
            continue
        rf = rows_f[s]
        cf = cols_f[s]
        r_mid = min(max(int(math.floor(rf)), 0), nrows - 1)
        c_mid = int(math.floor(cf))
        reach_r = int(math.ceil(radius_m[s] / dy_m))
        reach_c = int(math.ceil(radius_m[s] / dx_m[r_mid]))
        for r in range(max(r_mid - reach_r, 0), min(r_mid + reach_r + 1, nrows)):
            for c in range(max(c_mid - reach_c, 0), min(c_mid + reach_c + 1, ncols)):
                cell = r * ncols + c
                k = channel_id[cell]
                if k < 0 or excluded[cell]:
                    continue
                d = math.sqrt(((c + 0.5 - cf) * dx_m[r]) ** 2 + ((r + 0.5 - rf) * dy_m) ** 2)
                if d > radius_m[s]:
                    continue
                m = count[s]
                if m < k_max:
                    pos = m
                    count[s] = m + 1
                elif d < cand_dist[s, k_max - 1]:
                    pos = k_max - 1
                else:
                    continue
                while pos > 0 and cand_dist[s, pos - 1] > d:
                    cand[s, pos] = cand[s, pos - 1]
                    cand_dist[s, pos] = cand_dist[s, pos - 1]
                    pos -= 1
                cand[s, pos] = k
                cand_dist[s, pos] = d
    return cand, cand_dist, count


@njit(cache=True)
def _align_run(a: int, b: int, along: np.ndarray, cand: np.ndarray, cand_dist: np.ndarray, count: np.ndarray,
               tin: np.ndarray, tout: np.ndarray, to_outlet: np.ndarray, sigma: np.ndarray, beta: float,
               slack: float, max_gap: int, gap_penalty: float) -> tuple[float, int, int, int, int, int]:
    """
    The best local alignment of samples ``a..b-1`` onto one D8 path (Viterbi with free start and end).

    A state is a candidate channel cell of a sample. State ``k`` of sample ``i`` may follow state ``q`` of an
    earlier sample ``p`` (at most ``max_gap`` samples back) only when ``k`` is ``q`` or downstream of it, and then
    costs how far the D8 path from ``q`` to ``k`` differs in length from the source line between the two samples.
    A sample's fit is ``1 - (distance / sigma)^2 / 2``.

    Returns (score, first sample, last sample, first channel cell, last channel cell, matched samples).
    """
    n = b - a
    k_max = cand.shape[1]
    score = np.full((max(n, 1), k_max), -np.inf)
    back = np.full((max(n, 1), k_max), -1, np.int64)
    best = 0.0
    best_state = -1
    for ii in range(n):
        i = a + ii
        for j in range(count[i]):
            kj = cand[i, j]
            fit = 1.0 - 0.5 * (cand_dist[i, j] / sigma[i]) ** 2
            top = 0.0
            arg = -1
            for g in range(1, max_gap + 1):
                pp = ii - g
                if pp < 0:
                    break
                p = a + pp
                expected = along[i] - along[p]
                for q in range(count[p]):
                    sq = score[pp, q]
                    if sq <= 0.0:
                        continue
                    kq = cand[p, q]
                    if not (tin[kj] <= tin[kq] and tin[kq] <= tout[kj]):
                        continue
                    excess = abs(to_outlet[kq] - to_outlet[kj] - expected) - slack
                    if excess < 0.0:
                        excess = 0.0
                    s = sq - excess / beta - (g - 1) * gap_penalty
                    if s > top:
                        top = s
                        arg = pp * k_max + q
            score[ii, j] = top + fit
            back[ii, j] = arg
            if score[ii, j] > best:
                best = score[ii, j]
                best_state = ii * k_max + j
    if best_state < 0:
        return 0.0, -1, -1, -1, -1, 0
    last_sample = best_state // k_max
    last_cell = cand[a + last_sample, best_state % k_max]
    state = best_state
    matched = 0
    first_sample = last_sample
    first_cell = last_cell
    while state >= 0:
        ii = state // k_max
        matched += 1
        first_sample = ii
        first_cell = cand[a + ii, state % k_max]
        state = back[ii, state % k_max]
    return best, a + first_sample, a + last_sample, first_cell, last_cell, matched


@njit(cache=True)
def _align_reaches(sample_offsets: np.ndarray, usable: np.ndarray, along: np.ndarray, cand: np.ndarray,
                   cand_dist: np.ndarray, count: np.ndarray, tin: np.ndarray, tout: np.ndarray,
                   to_outlet: np.ndarray, sigma: np.ndarray, beta: float, slack: float, max_gap: int,
                   gap_penalty: float):
    """Align every reach: each run of usable samples separately, keeping the best-scoring run."""
    n_reaches = sample_offsets.size - 1
    score = np.zeros(n_reaches)
    first_sample = np.full(n_reaches, -1, np.int64)
    last_sample = np.full(n_reaches, -1, np.int64)
    first_cell = np.full(n_reaches, -1, np.int64)
    last_cell = np.full(n_reaches, -1, np.int64)
    matched = np.zeros(n_reaches, np.int64)
    for r in range(n_reaches):
        s = sample_offsets[r]
        stop = sample_offsets[r + 1]
        while s < stop:
            while s < stop and not usable[s]:
                s += 1
            run_start = s
            while s < stop and usable[s]:
                s += 1
            if s == run_start:
                continue
            res = _align_run(run_start, s, along, cand, cand_dist, count, tin, tout, to_outlet, sigma, beta,
                             slack, max_gap, gap_penalty)
            if res[0] > score[r]:
                score[r] = res[0]
                first_sample[r] = res[1]
                last_sample[r] = res[2]
                first_cell[r] = res[3]
                last_cell[r] = res[4]
                matched[r] = res[5]
    return score, first_sample, last_sample, first_cell, last_cell, matched


# ---------------------------------------------------------------------------------------------------------------
# Painting the matched reaches onto the D8 network
# ---------------------------------------------------------------------------------------------------------------

@njit(cache=True)
def _is_downstream(tin: np.ndarray, tout: np.ndarray, a: int, b: int) -> bool:
    """Whether node ``a`` is ``b`` or downstream of it, in a tour whose parents point downstream."""
    return tin[a] <= tin[b] and tin[b] <= tout[a]


@njit(cache=True)
def _related(r: int, x: int, ds: np.ndarray, stin: np.ndarray, stout: np.ndarray) -> bool:
    """
    Whether the source network sends reach ``r`` on to reach ``x``: x is downstream of r, or shares r's
    downstream reach (two tributaries of one river, or a tributary and a reach of that river).
    """
    if _is_downstream(stin, stout, x, r):
        return True
    return ds[r] >= 0 and ds[x] == ds[r]


@njit(cache=True)
def _extend_ends(start: np.ndarray, end: np.ndarray, down: np.ndarray, channel_id: np.ndarray,
                 to_outlet: np.ndarray, stop_mask: np.ndarray, ds: np.ndarray, stin: np.ndarray,
                 stout: np.ndarray, max_extension_m: float, snap_m: float) -> np.ndarray:
    """
    Where each matched reach's claim ends.

    A reach claims the D8 path from its matched start to its matched end, cut short at the first cell of
    ``stop_mask`` (a lake). Past its matched end it carries on down the D8 path, while that stays clear of
    ``stop_mask``, until it reaches a cell another reach matched onto: when the source network sends the reach on
    to that reach (see _related) it may go as far as ``max_extension_m`` to get there, and otherwise only
    ``snap_m``. A reach that reaches nobody ends where its match ended.
    """
    n = start.size
    # Which reaches matched onto each cell, as a sorted (cell, reach) list
    total = 0
    for r in range(n):
        if start[r] < 0:
            continue
        cell = start[r]
        while True:
            total += 1
            if cell == end[r] or down[cell] < 0:
                break
            cell = down[cell]
    owner_cell = np.empty(total, np.int64)
    owner = np.empty(total, np.int64)
    k = 0
    for r in range(n):
        if start[r] < 0:
            continue
        cell = start[r]
        while True:
            owner_cell[k] = cell
            owner[k] = r
            k += 1
            if cell == end[r] or down[cell] < 0:
                break
            cell = down[cell]
    order = np.argsort(owner_cell, kind="mergesort")
    owner_cell = owner_cell[order]
    owner = owner[order]

    claim_end = np.full(n, -1, np.int64)
    for r in range(n):
        if start[r] < 0:
            continue
        # The matched path, cut at the first stop cell
        cell = start[r]
        if stop_mask[cell]:
            continue
        last = cell
        while cell != end[r]:
            nxt = down[cell]
            if nxt < 0 or stop_mask[nxt]:
                break
            cell = nxt
            last = cell
        claim_end[r] = last
        if last != end[r]:
            continue
        # Past the matched end
        base = to_outlet[channel_id[last]]
        prev = last
        cell = down[last]
        while cell >= 0 and not stop_mask[cell]:
            travelled = base - to_outlet[channel_id[cell]]
            if travelled > max_extension_m:
                break
            lo = np.searchsorted(owner_cell, cell)
            hit = False
            joined = False
            while lo < total and owner_cell[lo] == cell:
                x = owner[lo]
                lo += 1
                if x == r:
                    continue
                hit = True
                if _related(r, x, ds, stin, stout) or travelled <= snap_m:
                    joined = True
            if hit:
                if joined:
                    claim_end[r] = prev
                break
            prev = cell
            cell = down[cell]
    return claim_end


@njit(cache=True)
def _outranks(a: int, b: int, order: np.ndarray, area: np.ndarray, score: np.ndarray) -> bool:
    if order[a] != order[b]:
        return order[a] > order[b]
    if area[a] != area[b]:
        return area[a] > area[b]
    if score[a] != score[b]:
        return score[a] > score[b]
    return a < b


@njit(cache=True)
def _paint(start: np.ndarray, claim_end: np.ndarray, down: np.ndarray, channel_id: np.ndarray,
           ctin: np.ndarray, ctout: np.ndarray, to_outlet: np.ndarray, ds: np.ndarray, stin: np.ndarray,
           stout: np.ndarray, order: np.ndarray, area: np.ndarray, score: np.ndarray, max_merge_m: float):
    """
    Label the claimed D8 paths with reaches, in flow order. See the module docstring for the rules.

    Returns the claimed cells (sorted), the reach on each (-1 for none), and each reach's first and last cell
    (-1 if it never won a cell).
    """
    n = start.size
    # The union of the claims. _extend_ends found each claim's end by walking this same path, so the walk
    # reaches it; the outlet check only guards against a corrupt pointer.
    total = 0
    for r in range(n):
        if claim_end[r] < 0:
            continue
        cell = start[r]
        while True:
            total += 1
            if cell == claim_end[r] or down[cell] < 0:
                break
            cell = down[cell]
    cells = np.empty(total, np.int64)
    k = 0
    for r in range(n):
        if claim_end[r] < 0:
            continue
        cell = start[r]
        while True:
            cells[k] = cell
            k += 1
            if cell == claim_end[r] or down[cell] < 0:
                break
            cell = down[cell]
    cells = np.unique(cells)
    m = cells.size

    # Upstream neighbours within the union, and the reaches starting at each cell
    downstream_index = np.full(m, -1, np.int64)
    pending = np.zeros(m, np.int64)
    for u in range(m):
        d = down[cells[u]]
        if d >= 0:
            j = np.searchsorted(cells, d)
            if j < m and cells[j] == d:
                downstream_index[u] = j
                pending[j] += 1
    first_up = np.zeros(m + 1, np.int64)
    for u in range(m):
        if downstream_index[u] >= 0:
            first_up[downstream_index[u] + 1] += 1
    for u in range(m):
        first_up[u + 1] += first_up[u]
    ups = np.empty(max(first_up[m], 1), np.int64)
    fill = first_up[:m].copy()
    for u in range(m):
        j = downstream_index[u]
        if j >= 0:
            ups[fill[j]] = u
            fill[j] += 1
    starters = np.empty(n, np.int64)
    starter_cell = np.empty(n, np.int64)
    n_starters = 0
    for r in range(n):
        if claim_end[r] >= 0:
            starters[n_starters] = r
            starter_cell[n_starters] = np.searchsorted(cells, start[r])
            n_starters += 1
    by_cell = np.argsort(starter_cell[:n_starters], kind="mergesort")
    starters = starters[:n_starters][by_cell]
    starter_cell = starter_cell[:n_starters][by_cell]

    label = np.full(m, -1, np.int64)
    first = np.full(n, -1, np.int64)
    last = np.full(n, -1, np.int64)
    alive = np.zeros(n, np.bool_)
    started = np.zeros(n, np.bool_)

    queue = np.empty(m, np.int64)
    head = 0
    tail = 0
    for u in range(m):
        if pending[u] == 0:
            queue[tail] = u
            tail += 1
    candidates = np.empty(n + 8, np.int64)
    while head < tail:
        u = queue[head]
        head += 1
        cell = cells[u]
        c = 0
        # Reaches flowing in from upstream that have not reached the end of their claim
        for t in range(first_up[u], first_up[u + 1]):
            v = ups[t]
            r = label[v]
            if r >= 0 and alive[r] and last[r] == cells[v] and cells[v] != claim_end[r]:
                candidates[c] = r
                c += 1
        # Reaches whose match starts here
        lo = np.searchsorted(starter_cell, u)
        while lo < n_starters and starter_cell[lo] == u:
            r = starters[lo]
            lo += 1
            if not started[r]:
                candidates[c] = r
                c += 1

        winner = -1
        if c == 1:
            winner = candidates[0]
        elif c > 1:
            # 1. A reach downstream of all the others in the source network carries on
            for i in range(c):
                x = candidates[i]
                ok = True
                for j in range(c):
                    if not _is_downstream(stin, stout, x, candidates[j]):
                        ok = False
                        break
                if ok:
                    winner = x
                    break
            # 2. Reaches that all drain into one source reach start it here, if its match lies further down
            #    this path
            if winner < 0:
                d = ds[candidates[0]]
                same = d >= 0
                for i in range(1, c):
                    if ds[candidates[i]] != d:
                        same = False
                        break
                if same and claim_end[d] >= 0 and not started[d] and channel_id[cell] >= 0 \
                        and channel_id[start[d]] >= 0:
                    cd = channel_id[start[d]]
                    cc = channel_id[cell]
                    if _is_downstream(ctin, ctout, cd, cc) and to_outlet[cc] - to_outlet[cd] <= max_merge_m:
                        # Every cell between here and its match must be claimed, or it would be cut in two
                        walk = cell
                        covered = True
                        while walk != start[d]:
                            walk = down[walk]
                            if walk < 0:
                                covered = False
                                break
                            j = np.searchsorted(cells, walk)
                            if j >= m or cells[j] != walk:
                                covered = False
                                break
                        if covered:
                            winner = d
            # 3. Otherwise the larger river carries on
            if winner < 0:
                winner = candidates[0]
                for i in range(1, c):
                    if _outranks(candidates[i], winner, order, area, score):
                        winner = candidates[i]
            for i in range(c):
                if candidates[i] != winner:
                    alive[candidates[i]] = False

        if winner >= 0:
            label[u] = winner
            if not started[winner]:
                started[winner] = True
                alive[winner] = True
                first[winner] = cell
            last[winner] = cell

        j = downstream_index[u]
        if j >= 0:
            pending[j] -= 1
            if pending[j] == 0:
                queue[tail] = j
                tail += 1
    return cells, label, first, last


# ---------------------------------------------------------------------------------------------------------------
# Putting it together
# ---------------------------------------------------------------------------------------------------------------

@dataclass
class ConflationResult:
    """The conflated network: one row per reach, and the stream raster it paints."""
    streams: gpd.GeoDataFrame
    raster: np.ndarray          # int32, the source id of each stream cell, 0 elsewhere
    matched_fraction: dict      # source id -> fraction of its usable length that was matched


def conflate_to_flow_directions(
        source_gdf: gpd.GeoDataFrame,
        flowdir: np.ndarray,
        geotransform: tuple,
        projection: str,
        id_col: str,
        ds_col: str,
        order_col: str | None = None,
        area_col: str | None = None,
        lake_mask: np.ndarray | None = None,
        valid_mask: np.ndarray | None = None,
        settings: ConflationSettings | None = None,
        flow: FlowNetwork | None = None,
) -> ConflationResult:
    """
    Move ``source_gdf``'s reaches onto the D8 network of ``flowdir`` (WhiteboxTools pointer codes).

    ``source_gdf`` must be in the raster's CRS. Its reaches keep their ids and every other column; each reach's
    line becomes the centres of the D8 cells it was matched onto, from upstream to downstream, and its downstream
    id becomes the reach its last cell drains into (-1 where that is no reach). Reaches that could not be matched
    are left out.

    ``lake_mask`` (True inside a lake) cuts the network at lakes: no reach is matched inside one, and a reach
    flowing into one ends at its shore. ``valid_mask`` (False on nodata or ocean) stops flow there.
    """
    settings = settings or ConflationSettings()
    if flow is None:
        flow = FlowNetwork(flowdir, geotransform, projection, valid_mask)
    nrows, ncols = flow.nrows, flow.ncols
    empty_raster = np.zeros((nrows, ncols), np.int32)
    out_columns = [c for c in source_gdf.columns if c != source_gdf.geometry.name]

    def empty_result() -> ConflationResult:
        return ConflationResult(gpd.GeoDataFrame(columns=out_columns + ["geometry"], geometry="geometry",
                                                 crs=source_gdf.crs), empty_raster, {})

    if source_gdf.empty:
        return empty_result()

    source = prepare_source_network(source_gdf, id_col, ds_col, order_col, area_col, flow)
    n = source.ids.size
    if n == 0:
        return empty_result()
    if source.ids.max() > np.iinfo(np.int32).max:
        raise ValueError(f"Stream ids up to {source.ids.max()} do not fit the stream raster's 32-bit cells.")

    stop = np.zeros(nrows * ncols, np.bool_) if lake_mask is None else \
        np.ascontiguousarray(lake_mask, dtype=np.bool_).ravel()

    # Sample each line in pixel space, scaled to metres so the spacing is even on a geographic grid
    gt = flow.geotransform
    coords, index = shapely.get_coordinates(source.geometry, return_index=True)
    cols_f = (coords[:, 0] - gt[0]) / gt[1]
    rows_f = (coords[:, 1] - gt[3]) / gt[5]
    dx_mid = float(np.median(flow.dx_m))
    offsets = np.searchsorted(index, np.arange(n + 1)).astype(np.int64)
    sx, sy, along, sample_offsets = _resample(cols_f * dx_mid, rows_f * flow.dy_m, offsets,
                                              settings.sample_spacing_cells * flow.cell_size_m)
    s_cols = sx / dx_mid
    s_rows = sy / flow.dy_m
    reach_of_sample = np.repeat(np.arange(n), np.diff(sample_offsets))

    inside = (s_rows >= 0) & (s_rows < nrows) & (s_cols >= 0) & (s_cols < ncols)
    cell_of_sample = np.where(inside, np.floor(s_rows).clip(0, nrows - 1).astype(np.int64) * ncols
                              + np.floor(s_cols).clip(0, ncols - 1).astype(np.int64), 0)
    usable = inside & flow.valid.ravel()[cell_of_sample] & ~stop[cell_of_sample]

    flow.build_channels(settings.min_channel_area_km2, _inflow_cells(sample_offsets, usable, cell_of_sample))
    if flow.channel_cells.size == 0:
        return empty_result()

    # Search wider, and fit looser, along big rivers
    reach_area = np.where(source.area > 0, source.area, 0.0)
    reach_area = np.maximum(reach_area, _approximate_areas(source, flow))
    radius = np.clip(settings.search_radius_cells * flow.cell_size_m
                     + settings.search_radius_per_sqrt_km2 * np.sqrt(reach_area),
                     settings.search_radius_cells * flow.cell_size_m, settings.max_search_radius_m)
    sigma = settings.sigma_cells * flow.cell_size_m * radius / (settings.search_radius_cells * flow.cell_size_m)
    sample_radius = radius[reach_of_sample]
    sample_sigma = sigma[reach_of_sample]

    cand, cand_dist, count = _find_candidates(s_rows, s_cols, usable, sample_radius, settings.max_candidates,
                                              flow.channel_id, stop, nrows, ncols, flow.dx_m, flow.dy_m)
    # How much of each reach could have been matched at all: inside the domain and outside lakes
    usable_length = np.bincount(reach_of_sample[usable], minlength=n).astype(np.float64)
    usable &= count > 0
    score, first_sample, last_sample, first_cid, last_cid, matched = _align_reaches(
        sample_offsets, usable, along, cand, cand_dist, count, flow.tin, flow.tout, flow.to_outlet_m, sample_sigma,
        settings.beta_cells * flow.cell_size_m, settings.slack_cells * flow.cell_size_m, settings.max_gap_samples,
        settings.gap_penalty)

    matched_fraction = np.where((matched > 0) & (usable_length > 0),
                                (last_sample - first_sample + 1) / np.maximum(usable_length, 1), 0.0)
    keep = (matched >= settings.min_matched_samples) & (matched_fraction >= settings.min_matched_fraction)
    matched_fraction[~keep] = 0.0
    start = np.where(keep, flow.channel_cells[np.maximum(first_cid, 0)], -1)
    end = np.where(keep, flow.channel_cells[np.maximum(last_cid, 0)], -1)

    claim_end = _extend_ends(start, end, flow.down, flow.channel_id, flow.to_outlet_m, stop, source.ds,
                             source.tin, source.tout, settings.max_extension_m,
                             settings.snap_cells * flow.cell_size_m)
    cells, label, first, last = _paint(start, claim_end, flow.down, flow.channel_id, flow.tin, flow.tout,
                                       flow.to_outlet_m, source.ds, source.tin, source.tout, source.order,
                                       reach_area, score, settings.max_extension_m)
    label = _fold_short_reaches(cells, label, first, last, flow, settings.min_reach_cells, source, reach_area,
                                score)

    raster = np.zeros(nrows * ncols, np.int32)
    labelled = label >= 0
    raster[cells[labelled]] = source.ids[label[labelled]]
    raster = raster.reshape(nrows, ncols)

    streams = _reach_lines(cells, label, flow, source, source_gdf, id_col, ds_col)
    LOG.info(
        f"Conflation: {len(streams)} of {n} source reaches placed on the D8 network "
        f"({int(keep.sum())} matched, median matched fraction {np.median(matched_fraction[keep]) if keep.any() else 0:.2f})."
    )
    fraction = {int(source.ids[i]): float(matched_fraction[i]) for i in range(n)}
    return ConflationResult(streams, raster, fraction)


def _inflow_cells(sample_offsets: np.ndarray, usable: np.ndarray, cell_of_sample: np.ndarray) -> np.ndarray:
    """
    Where water the domain's DEM doesn't drain flows in along a source reach: the first usable sample of each
    reach that starts outside the domain, on nodata or in a lake.
    """
    cells = []
    for r in range(sample_offsets.size - 1):
        a, b = sample_offsets[r], sample_offsets[r + 1]
        ok = np.flatnonzero(usable[a:b])
        if ok.size and ok[0] > 0:
            cells.append(cell_of_sample[a + ok[0]])
    return np.unique(np.asarray(cells, np.int64))


def _approximate_areas(source: SourceNetwork, flow: FlowNetwork) -> np.ndarray:
    """Drainage area at each reach's downstream end, read off the DEM when the source network has none."""
    if (source.area > 0).all():
        return np.zeros(source.ids.size)
    ends = shapely.points([g.coords[-1] for g in source.geometry])
    return np.array([_area_near(flow, p) for p in ends], np.float64)


def _fold_short_reaches(cells: np.ndarray, label: np.ndarray, first: np.ndarray, last: np.ndarray,
                        flow: FlowNetwork, min_cells: int, source: SourceNetwork, area: np.ndarray,
                        score: np.ndarray) -> np.ndarray:
    """
    Give the cells of every reach shorter than ``min_cells`` to the reach that flows into its first cell (the
    largest, if several do), or unlabel them when none does. A reach is at least two cells, the fewest a line can
    be drawn through.
    """
    min_cells = max(min_cells, 2)
    counts = np.bincount(label[label >= 0], minlength=first.size)
    short = np.flatnonzero((counts > 0) & (counts < min_cells))
    if short.size == 0:
        return label
    down = flow.down
    label = label.copy()
    index = {int(c): i for i, c in enumerate(cells)}
    # Upstream first, so that a short reach below another short reach finds it already folded into its feeder
    short = sorted(short, key=lambda r: -flow.to_outlet_m[flow.channel_id[first[r]]])
    for r in short:
        head = int(first[r])
        feeders = [int(x) for x in np.flatnonzero(last >= 0)
                   if x != r and counts[x] >= min_cells and down[last[x]] == head]
        heir = -1
        for x in feeders:
            if heir < 0 or (source.order[x], area[x], score[x]) > (source.order[heir], area[heir], score[heir]):
                heir = x
        cell = head
        while True:
            label[index[cell]] = heir
            if cell == last[r]:
                break
            cell = int(down[cell])
        if heir >= 0:
            last[heir] = last[r]
            counts[heir] += counts[r]
        first[r] = last[r] = -1
        counts[r] = 0
    return label


def _reach_lines(cells: np.ndarray, label: np.ndarray, flow: FlowNetwork, source: SourceNetwork,
                 source_gdf: gpd.GeoDataFrame, id_col: str, ds_col: str) -> gpd.GeoDataFrame:
    """One row per painted reach: its cells as a line, upstream to downstream, and the reach it drains into."""
    labelled = label >= 0
    if not labelled.any():
        columns = [c for c in source_gdf.columns if c != source_gdf.geometry.name]
        return gpd.GeoDataFrame(columns=columns + ["geometry"], geometry="geometry", crs=source_gdf.crs)
    reach_cells = cells[labelled]
    reach_label = label[labelled]
    # Order each reach's cells from upstream to downstream: along a D8 path the distance to the outlet falls
    to_outlet = flow.to_outlet_m[flow.channel_id[reach_cells]]
    order = np.lexsort((-to_outlet, reach_label))
    reach_cells = reach_cells[order]
    reach_label = reach_label[order]
    bounds = np.flatnonzero(np.diff(reach_label)) + 1
    starts = np.concatenate(([0], bounds))
    stops = np.concatenate((bounds, [reach_label.size]))

    cell_label = dict(zip(cells[labelled].tolist(), label[labelled].tolist()))
    xy = flow.cell_centres(reach_cells)
    rows = []
    for a, b in zip(starts, stops):
        r = int(reach_label[a])
        downstream = int(flow.down[reach_cells[b - 1]])
        ds_reach = cell_label.get(downstream, -1) if downstream >= 0 else -1
        ds_id = int(source.ids[ds_reach]) if ds_reach >= 0 and ds_reach != r else -1
        rows.append((int(source.ids[r]), ds_id, shapely.linestrings(xy[a:b])))
    out = pd.DataFrame(rows, columns=[id_col, "__ds", "geometry"])
    out["topological_order"] = _upstream_first(out[id_col].to_numpy(), out["__ds"].to_numpy())

    attributes = source_gdf.drop(columns=[source_gdf.geometry.name]).copy()
    attributes = attributes.drop(columns=[c for c in ("topological_order",) if c in attributes.columns])
    attributes[id_col] = pd.to_numeric(attributes[id_col], errors="coerce")
    attributes = attributes.drop_duplicates(subset=id_col).set_index(id_col)
    out = out.join(attributes, on=id_col)
    out[ds_col] = out.pop("__ds")
    columns = [id_col, ds_col] + [c for c in out.columns if c not in (id_col, ds_col, "geometry")] + ["geometry"]
    return gpd.GeoDataFrame(out[columns], geometry="geometry", crs=source_gdf.crs)


def _upstream_first(ids: np.ndarray, ds: np.ndarray) -> np.ndarray:
    """
    A rank for each reach that puts every reach before the one it drains into; ties go to the smaller id.

    Curve2Flood's FLDPLN mapper sorts the stream network by this ``topological_order`` before chaining reaches
    into main stems.
    """
    index = {int(rid): i for i, rid in enumerate(ids)}
    downstream = np.array([index.get(int(d), -1) for d in ds], np.int64)
    pending = np.zeros(ids.size, np.int64)
    for d in downstream:
        if d >= 0:
            pending[d] += 1
    ready = [(int(ids[i]), i) for i in range(ids.size) if pending[i] == 0]
    heapq.heapify(ready)
    rank = np.zeros(ids.size, np.int64)
    k = 0
    while ready:
        _, i = heapq.heappop(ready)
        rank[i] = k
        k += 1
        d = downstream[i]
        if d >= 0:
            pending[d] -= 1
            if pending[d] == 0:
                heapq.heappush(ready, (int(ids[d]), d))
    return rank
