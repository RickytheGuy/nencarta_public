from pathlib import Path

import pytest

from nencarta.core.configs import NencartaConfig
from nencarta.core.defaults import DEFAULT_CONFIG
from nencarta.core.enumerations import FloodMapMode, Mapper, StreamflowSource
from nencarta.workspace import Workspace


def base_config(workspace_tmp: Path, **overrides) -> dict:
    dem_dir = workspace_tmp / "dems"
    dem_dir.mkdir(exist_ok=True)
    config = {
        "name": "example",
        "flowline": str(workspace_tmp / "streams.parquet"),
        "dem_dir": str(dem_dir),
        "output_dir": str(workspace_tmp / "out"),
        "mapper": "Curve2Flood-Kernel Weighted",
        "streamflow_source": "GEOGLOWS",
        "use_specified_depth_for_bathy_mask": True,
        "specify_depths_for_bathy_mask": [0.1],
        "clean_dem": False,
        "make_fist_inputs": False,
    }
    config.update(overrides)
    return config


def test_config_normalizes_core_options(workspace_tmp):
    config = NencartaConfig(base_config(workspace_tmp))

    assert config.watershed_name == "example"
    assert config.mapper == Mapper.CURVE2FLOOD_KERNEL_WEIGHTED
    assert config.streamflow_source == StreamflowSource.GEOGLOWS
    assert config.floodmap_mode == FloodMapMode.FORECAST
    assert config.dem_filter == "*"
    assert config.overwrite is False


def test_config_accepts_overwrite_option(workspace_tmp):
    config = NencartaConfig(base_config(workspace_tmp, overwrite=True))

    assert config.overwrite is True


def test_config_uses_shared_defaults_for_bathymetry_options(workspace_tmp):
    config = NencartaConfig(base_config(workspace_tmp))

    assert config.use_dem_derived_channel_mask is DEFAULT_CONFIG["use_dem_derived_channel_mask"]
    assert config.coefficient_depth == DEFAULT_CONFIG["coefficient_depth"]
    assert config.exponent_depth == DEFAULT_CONFIG["exponent_depth"]
    assert config.coefficient_width == DEFAULT_CONFIG["coefficient_width"]
    assert config.exponent_width == DEFAULT_CONFIG["exponent_width"]


def test_clean_dem_requires_two_bathy_mask_depths(workspace_tmp):
    with pytest.raises(ValueError, match="specify_depths_for_bathy_mask"):
        NencartaConfig(base_config(workspace_tmp, clean_dem=True))


def test_workspace_uses_expected_output_paths(workspace_tmp):
    dem = workspace_tmp / "dems" / "fabdem.tif"
    dem.parent.mkdir(exist_ok=True)
    dem.touch()

    workspace = Workspace(NencartaConfig(base_config(workspace_tmp)), dem)

    assert workspace.output_dir == workspace_tmp / "out" / "example"
    assert workspace.FileName == "fabdem"
    assert workspace.assigned_dem == dem
    assert workspace.DEM_StrmShp.name == "GEOGLOWS_fabdem_StrmShp.gpkg"
    assert workspace.FloodWSEFile.name == "GEOGLOWS_fabdem_ARC_FloodWSE.tif"


def test_workspace_bbox_names_synthetic_dem(workspace_tmp):
    config = base_config(
        workspace_tmp,
        dem_dir=None,
        source_dems=[str(workspace_tmp / "source.tif")],
        bbox=[-105.123456, 40.123456, -104.123456, 41.123456],
        buffer=True,
    )

    workspace = Workspace(NencartaConfig(config), None)

    assert workspace.FileName == "dem_-105.12346_40.12346_-104.12346_41.12346_buffered"
    assert workspace.assigned_dem.name.endswith(".tif")


def test_short_file_names_accept_a_dem_used_as_is(workspace_tmp):
    # Without a buffer or bbox the assigned DEM is the input DEM itself, so the two share a path by design
    dem = workspace_tmp / "dems" / "fabdem.tif"
    dem.parent.mkdir(exist_ok=True)
    dem.touch()

    workspace = Workspace(NencartaConfig(base_config(workspace_tmp, short_file_names=True)), dem)

    assert workspace.original_dem == workspace.assigned_dem == dem


def test_file_names_that_give_two_files_one_path_raise(workspace_tmp):
    dem = workspace_tmp / "dems" / "fabdem.tif"
    dem.parent.mkdir(exist_ok=True)
    dem.touch()
    config = base_config(workspace_tmp, short_file_names=True, file_names={"fixed": "dem_out", "filled": "dem_out"})

    with pytest.raises(ValueError, match="same path"):
        Workspace(NencartaConfig(config), dem)


def test_a_file_named_onto_the_input_dem_raises(workspace_tmp):
    # The DEM may share its path with the assigned DEM, but an output written over it would destroy the input
    dem = workspace_tmp / "dems" / "fabdem.tif"
    dem.parent.mkdir(exist_ok=True)
    dem.touch()
    config = base_config(workspace_tmp, short_file_names=True, file_names={"fixed": "fabdem"},
                         folder_paths={"DEM_Updated": str(workspace_tmp / "dems")})

    with pytest.raises(ValueError, match="same path"):
        Workspace(NencartaConfig(config), dem)
