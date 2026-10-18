"""
Where a DEM shows its rivers' channels, whatever the DEM.

Some DEMs flatten their water (FABDEM and Copernicus to whole half metres, USGS 3DEP to any height), others store
heights in whole metres (SRTM, NASADEM, ASTER, AW3D30) or leave the water as noisy as the land. Two signs of a channel
hold across them, given the streams:

- Flat water: cells exactly as high as at least two of their neighbours, in patches of more than a few cells, no
  more than ``flat_rise_m`` above the river beside them, on a DEM whose heights are rounded to a centimetre or finer
  (FABDEM, Copernicus, USGS 3DEP). On a DEM rounded more coarsely (GEDTM's 0.1 m, SRTM's whole metres) most of the
  land is flat too, and flatness says nothing.
- Low ground by the stream: no more than ``rise_m`` above the river beside it, within ``widths`` of the river's
  expected width (``width_per_sqrt_km2`` metres per square root of the drainage area in km2) of the stream, along a
  river expected to be at least a cell wide. A narrower river's channel is its stream cells.

The river beside a cell is its nearest stream cell, whose level is the lowest of the DEM within ``snap_m`` of it (the
stream may be drawn a cell or two off the water), averaged along the streams over ``smooth_m``. On the eight DEMs of the
FIM sites the two together find ESA's river water along big rivers with a median F1 of 0.63 to 0.82 on six (0.33 on
AW3D30, 0.20 on GEDTM), against 0.84 on FABDEM and Copernicus, nothing on USGS 3DEP and about 0.1 elsewhere for whole
half metres alone.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numba import njit
from scipy.ndimage import binary_dilation, distance_transform_edt, label

# Rows of the grid worked on at a time where every cell's nearest stream cell is used, so the temporary arrays stay small
_CHUNK_ROWS = 512


@dataclass
class ChannelSettings:
    # How far from a stream a channel is looked for at all
    corridor_m: float = 1500.0
    # The river's level: the lowest of the DEM within snap_m of a stream cell, averaged along the streams over smooth_m
    snap_m: float = 60.0
    smooth_m: float = 150.0
    # Flat water: a cell as high as at least two neighbours, in a patch of more than flat_min_cells such cells, no more
    # than flat_rise_m above the river, on a DEM whose heights are rounded to no more than max_flat_step_m (see
    # height_step)
    flat_min_cells: int = 6
    flat_rise_m: float = 2.0
    max_flat_step_m: float = 0.01
    # Low ground: no more than rise_m above the river, within widths times the river's expected width of the stream,
    # along rivers expected to be at least min_width_cells wide (a narrower one's channel is its stream cells)
    rise_m: float = 1.0
    widths: float = 1.5
    width_per_sqrt_km2: float = 1.1
    min_width_cells: float = 1.0


@dataclass
class Channel:
    """
    The flat water and the whole channel (flat water or low ground by the streams), as boolean rasters, or None for
    one river_channel was told not to work out.
    """
    flat_water: np.ndarray | None
    channel: np.ndarray | None


def flat_cells(dem: np.ndarray, valid: np.ndarray, min_cells: int = 6) -> np.ndarray:
    """Cells exactly as high as at least two of their neighbours, in 8-connected patches of more than ``min_cells``."""
    flat = _level_with_two_neighbours(np.asarray(dem), np.ascontiguousarray(valid, dtype=np.bool_))
    labels, _ = label(flat, np.ones((3, 3), int))
    # Only the few flat cells are counted and looked up, not the whole grid
    patch = labels[flat]
    out = np.zeros(flat.shape, bool)
    out[flat] = np.bincount(patch)[patch] > min_cells
    return out


@njit(cache=True)
def _level_with_two_neighbours(dem: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """The valid cells exactly as high as at least two of their valid neighbours."""
    nrows, ncols = dem.shape
    out = np.zeros((nrows, ncols), np.bool_)
    for r0 in range(nrows):
        for c0 in range(ncols):
            if not valid[r0, c0]:
                continue
            # The cell itself counts too, so three in its 3 x 3 block (none if it is NaN, which equals nothing)
            z = dem[r0, c0]
            equal = 0
            if 0 < r0 < nrows - 1 and 0 < c0 < ncols - 1:
                for r in (r0 - 1, r0, r0 + 1):
                    for c in (c0 - 1, c0, c0 + 1):
                        if dem[r, c] == z and valid[r, c]:
                            equal += 1
            else:
                for r in range(max(r0 - 1, 0), min(r0 + 2, nrows)):
                    for c in range(max(c0 - 1, 0), min(c0 + 2, ncols)):
                        if dem[r, c] == z and valid[r, c]:
                            equal += 1
            out[r0, c0] = equal >= 3
    return out


def height_step(dem: np.ndarray, valid: np.ndarray, sample: int = 200_000) -> float:
    """
    The step a DEM's heights are rounded to: the coarsest of 1 m, 0.1 m, 1 cm and 1 mm that at least 90% of its cells
    (up to ``sample`` of them) are whole multiples of, or 0 when they are continuous.
    """
    z = np.asarray(dem)[valid].astype(np.float32)
    if z.size > sample:
        z = z[np.random.default_rng(0).choice(z.size, sample, replace=False)]
    if z.size == 0:
        return 0.0
    for step in (1.0, 0.1, 0.01, 0.001):
        if np.mean(np.float32(np.round(z.astype(np.float64) / step) * step) == z) >= 0.9:
            return step
    return 0.0


def touching(mask: np.ndarray, streams: np.ndarray) -> np.ndarray:
    """The 8-connected patches of ``mask`` that touch a stream cell, or a cell next to one."""
    near = binary_dilation(streams > 0, np.ones((3, 3), bool))
    labels, _ = label(mask, np.ones((3, 3), int))
    keep = np.zeros(labels.max() + 1, bool)
    keep[np.unique(labels[near & (labels > 0)])] = True
    keep[0] = False
    return keep[labels]


def river_channel(dem: np.ndarray, valid: np.ndarray, streams: np.ndarray, area_km2: dict, dx_m: float, dy_m: float,
                  settings: ChannelSettings | None = None, connected: bool = True, flat_water: bool = True,
                  channel: bool = True, double: bool = False) -> Channel:
    """
    The channel the DEM shows along ``streams`` (each stream cell's reach id, 0 elsewhere): see the module. ``area_km2``
    is each reach's drainage area, by id (a reach without one is taken as draining 1 km2); ``dx_m`` and ``dy_m`` are
    the cell's size. With ``connected``, only the patches touching a stream are kept. With ``flat_water`` or
    ``channel`` False, that one isn't worked out, and is None. With ``double``, the river's level is averaged in double
    precision on a single-precision DEM too, as on a float64 copy of it (see _stream_levels).
    """
    s = settings or ChannelSettings()
    valid = np.ascontiguousarray(valid, dtype=np.bool_)
    streams = np.asarray(streams)
    on_stream = valid & (streams > 0)
    empty = np.zeros(dem.shape, bool)
    if not on_stream.any():
        return Channel(empty if flat_water else None, empty.copy() if channel else None)
    cell = math.sqrt(dx_m * dy_m)
    flat = flat_cells(dem, valid, s.flat_min_cells) if height_step(dem, valid) <= s.max_flat_step_m else None
    if not channel and flat is None:
        return Channel(empty, None)

    # The river's level at each stream cell, and the river beside each cell: its nearest stream cell
    rows, cols, level = _stream_levels(dem, valid, on_stream, cell, s, double)
    have = np.zeros(dem.shape, bool)
    have[rows, cols] = True
    sampling = (dy_m, dx_m)
    nearest = distance_transform_edt(~have, sampling=sampling, return_distances=False, return_indices=True)
    del have
    index = np.full(dem.shape, -1, np.int32)
    index[rows, cols] = np.arange(rows.size, dtype=np.int32)

    if not channel:
        # Only the flat cells' rise above the river is needed
        flat_rows, flat_cols = np.nonzero(flat)
        near_rows, near_cols = nearest[0][flat_rows, flat_cols], nearest[1][flat_rows, flat_cols]
        del nearest
        distance = _feature_distance(near_rows - flat_rows, near_cols - flat_cols, sampling)
        rise = dem[flat_rows, flat_cols] - level[index[near_rows, near_cols]]
        keep = (distance <= s.corridor_m) & (rise <= s.flat_rise_m)
        water = np.zeros(dem.shape, bool)
        water[flat_rows[keep], flat_cols[keep]] = True
        return Channel(touching(water, on_stream) if connected else water, None)

    ids = np.unique(streams[rows, cols])
    width = s.width_per_sqrt_km2 * np.sqrt([max(area_km2.get(int(i), 1.0), 0.01) for i in ids])
    width = np.where(width >= s.min_width_cells * cell, width, 0.0)[np.searchsorted(ids, streams[rows, cols])]
    water = np.zeros(dem.shape, bool)
    whole = np.zeros(dem.shape, bool)
    ncols = dem.shape[1]
    for r0 in range(0, dem.shape[0], _CHUNK_ROWS):
        r1 = min(r0 + _CHUNK_ROWS, dem.shape[0])
        near_rows, near_cols = nearest[0, r0:r1], nearest[1, r0:r1]
        distance = _feature_distance(near_rows - np.arange(r0, r1, dtype=np.int32)[:, None],
                                     near_cols - np.arange(ncols, dtype=np.int32)[None, :], sampling)
        k = index[near_rows, near_cols]
        rise = dem[r0:r1] - level[k]
        corridor = valid[r0:r1] & (distance <= s.corridor_m)
        if flat is not None:
            water[r0:r1] = flat[r0:r1] & corridor & (rise <= s.flat_rise_m)
        w = width[k]
        whole[r0:r1] = water[r0:r1] | (corridor & (rise <= s.rise_m) & (w > 0) & (distance <= s.widths * w))
    del nearest, index
    if connected:
        water = touching(water, on_stream) if flat_water else None
        whole = touching(whole, on_stream)
    return Channel(water if flat_water else None, whole)


def _feature_distance(d_row: np.ndarray, d_col: np.ndarray, sampling: tuple[float, float]) -> np.ndarray:
    """
    The distance in metres to a cell ``d_row`` rows and ``d_col`` columns away (``sampling``: the cells' height and
    width), worked out as scipy's distance_transform_edt works it out from its feature transform, to the last bit.
    """
    d_row = d_row.astype(np.float64)
    d_col = d_col.astype(np.float64)
    d_row *= sampling[0]
    d_col *= sampling[1]
    d_row *= d_row
    d_col *= d_col
    d_row += d_col
    return np.sqrt(d_row, out=d_row)


def _stream_levels(dem: np.ndarray, valid: np.ndarray, on_stream: np.ndarray, cell: float, s: ChannelSettings,
                   double: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    The river's level at the stream cells: the lowest of the DEM within snap_m of each, averaged along the streams by a
    Gaussian of smooth_m. Returns the stream cells' rows and columns, and their levels.

    Worked out at the stream cells alone, but exactly as a minimum filter and scipy's gaussian_filter (reflected at the
    edges, truncated at 4 sigma) over the whole grid work it out, to the last bit: in the DEM's precision when that is
    single (see _smoothed_at), unless ``double``.
    """
    k = max(int(round(s.snap_m / cell)), 1)
    rows, cols = np.nonzero(on_stream)
    lowest = _window_min(dem, valid, rows, cols, k)
    finite = np.isfinite(lowest)
    rows, cols, lowest = rows[finite], cols[finite], lowest[finite]

    # gaussian_filter's kernel: truncated at 4 sigma, normalised to sum to 1
    sigma = float(max(s.smooth_m / cell, 1.0))
    radius = int(4.0 * sigma + 0.5)
    x = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 / (sigma * sigma) * x ** 2)
    kernel = np.ascontiguousarray((kernel / kernel.sum())[::-1])
    # A minimum filter over the DEM works in its type, or in float64 if that holds no infinity. The lowest heights are
    # the DEM's own, so a single-precision DEM's are held in single precision whichever the filters work in.
    work = np.where(np.ones(1, bool), np.zeros(1, dem.dtype), np.inf).dtype
    values = np.zeros(dem.shape, np.float32 if dem.dtype == np.float32 else np.float64)
    values[rows, cols] = lowest
    weight = np.zeros(dem.shape, np.bool_)
    weight[rows, cols] = True
    return rows, cols, _smoothed_at(values, weight, rows, cols, kernel, work == np.float32 and not double)


@njit(cache=True)
def _window_min(dem: np.ndarray, valid: np.ndarray, rows: np.ndarray, cols: np.ndarray, k: int) -> np.ndarray:
    """The lowest valid cell of ``dem`` within ``k`` cells (a square) of each of the cells, inf where there is none."""
    nrows, ncols = dem.shape
    out = np.full(rows.size, np.inf)
    for i in range(rows.size):
        for r in range(max(rows[i] - k, 0), min(rows[i] + k + 1, nrows)):
            for c in range(max(cols[i] - k, 0), min(cols[i] + k + 1, ncols)):
                if valid[r, c] and dem[r, c] < out[i]:
                    out[i] = dem[r, c]
    return out


@njit(cache=True)
def _reflect(i: int, n: int) -> int:
    """Index ``i`` of a line of ``n`` reflected back into it about its ends (d c b a | a b c d | d c b a)."""
    m = i % (2 * n)
    return m if m < n else 2 * n - 1 - m


@njit(cache=True)
def _smoothed_at(values: np.ndarray, weight: np.ndarray, rows: np.ndarray, cols: np.ndarray, kernel: np.ndarray,
                 single: bool) -> np.ndarray:
    """
    At each of the cells, gaussian_filter(values) / max(gaussian_filter(weight), 1e-12), with scipy's symmetric
    ``kernel`` and its reflected edges, summed in the order scipy sums them (down the columns, then along the rows,
    each pair of cells either side of the middle from the outermost in) so the result is the same to the last bit.
    With ``single``, ``values`` are float32 and each pass is rounded to float32, as scipy stores them.
    """
    nrows, ncols = values.shape
    half = kernel.size // 2
    middle = kernel[half]
    down_v = np.empty(kernel.size)
    down_w = np.empty(kernel.size)
    out = np.empty(rows.size)
    for i in range(rows.size):
        r, c = rows[i], cols[i]
        inside_rows = half <= r < nrows - half
        for j in range(-half, half + 1):
            col = c + j if 0 <= c + j < ncols else _reflect(c + j, ncols)
            v = np.float64(values[r, col]) * middle
            w = (1.0 if weight[r, col] else 0.0) * middle
            for jj in range(-half, 0):
                if inside_rows:
                    above, below = r + jj, r - jj
                else:
                    above, below = _reflect(r + jj, nrows), _reflect(r - jj, nrows)
                v += (np.float64(values[above, col]) + np.float64(values[below, col])) * kernel[half + jj]
                w += ((1.0 if weight[above, col] else 0.0) + (1.0 if weight[below, col] else 0.0)) * kernel[half + jj]
            down_v[j + half] = np.float64(np.float32(v)) if single else v
            down_w[j + half] = w
        total = down_v[half] * middle
        norm = down_w[half] * middle
        for jj in range(-half, 0):
            total += (down_v[half + jj] + down_v[half - jj]) * kernel[half + jj]
            norm += (down_w[half + jj] + down_w[half - jj]) * kernel[half + jj]
        if single:
            total = np.float64(np.float32(total))
        out[i] = total / max(norm, 1e-12)
    return out
