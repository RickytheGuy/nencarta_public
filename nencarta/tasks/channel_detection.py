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
from scipy.ndimage import binary_dilation, distance_transform_edt, gaussian_filter, label, minimum_filter


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
    """The flat water and the whole channel (flat water or low ground by the streams), as boolean rasters."""
    flat_water: np.ndarray
    channel: np.ndarray


def flat_cells(dem: np.ndarray, valid: np.ndarray, min_cells: int = 6) -> np.ndarray:
    """Cells exactly as high as at least two of their neighbours, in 8-connected patches of more than ``min_cells``."""
    z = np.where(valid, dem, np.nan)
    padded = np.pad(z, 1, constant_values=np.nan)
    equal = np.zeros(z.shape, np.int8)
    for d_row in (-1, 0, 1):
        for d_col in (-1, 0, 1):
            if d_row or d_col:
                equal += padded[1 + d_row:1 + d_row + z.shape[0], 1 + d_col:1 + d_col + z.shape[1]] == z
    flat = valid & (equal >= 2)
    labels, _ = label(flat, np.ones((3, 3), int))
    size = np.bincount(labels.ravel())
    size[0] = 0
    return size[labels] > min_cells


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
                  settings: ChannelSettings | None = None, connected: bool = True) -> Channel:
    """
    The channel the DEM shows along ``streams`` (each stream cell's reach id, 0 elsewhere): see the module. ``area_km2``
    is each reach's drainage area, by id (a reach without one is taken as draining 1 km2); ``dx_m`` and ``dy_m`` are
    the cell's size. With ``connected``, only the patches touching a stream are kept.
    """
    s = settings or ChannelSettings()
    valid = np.asarray(valid, bool)
    streams = np.where(valid, np.asarray(streams), 0)
    on_stream = streams > 0
    empty = np.zeros(dem.shape, bool)
    if not on_stream.any():
        return Channel(empty, empty.copy())
    cell = math.sqrt(dx_m * dy_m)

    # The river's level at each stream cell
    k = max(int(round(s.snap_m / cell)), 1)
    lowest = minimum_filter(np.where(valid, dem, np.inf), size=2 * k + 1)
    have = on_stream & np.isfinite(lowest)
    weight = have.astype(np.float64)
    sigma = max(s.smooth_m / cell, 1.0)
    total = gaussian_filter(np.where(have, lowest, 0.0), sigma)
    norm = gaussian_filter(weight, sigma)
    level = np.full(dem.shape, np.nan)
    level[have] = total[have] / np.maximum(norm[have], 1e-12)
    del lowest, total, norm, weight

    # The river beside each cell: its nearest stream cell
    distance, (rows, cols) = distance_transform_edt(~have, sampling=(dy_m, dx_m), return_indices=True)
    rise = dem - level[rows, cols]
    owner = streams[rows, cols]
    del rows, cols
    corridor = valid & (distance <= s.corridor_m)
    rise = np.where(corridor, rise, np.inf)

    ids, inverse = np.unique(owner, return_inverse=True)
    width = s.width_per_sqrt_km2 * np.sqrt([max(area_km2.get(int(i), 1.0), 0.01) for i in ids])
    width = np.where(width >= s.min_width_cells * cell, width, 0.0)[inverse.reshape(owner.shape)]
    near = (width > 0) & (distance <= s.widths * width)
    del owner, inverse, distance

    if height_step(dem, valid) <= s.max_flat_step_m:
        flat_water = flat_cells(dem, valid, s.flat_min_cells) & (rise <= s.flat_rise_m)
    else:
        flat_water = empty.copy()
    channel = flat_water | ((rise <= s.rise_m) & near)
    if connected:
        flat_water = touching(flat_water, streams)
        channel = touching(channel, streams)
    return Channel(flat_water, channel)
