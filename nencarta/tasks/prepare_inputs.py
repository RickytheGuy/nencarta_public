
from nencarta.tasks import (
    assign_and_validate_dem,
    Create_BaseLine_Manning_n_File_ESA,
    make_land_cover,
    make_stream_geometry,
    make_fldpln_inputs,
    define_configs_for_dem_cleaning,
    make_clean_dem,
    make_reanalysis_file,
    make_stream_raster,
    make_flood_flow_file_from_base_max_file,
    make_flow_file_from_forecast,
    assign_user_flow_files,
    make_return_period_flow_file,
    define_arc_configs,
    define_mapper_configs,
    make_water_mask,
    remove_old_forecast_files,
    get_fist_inputs)
from nencarta.workspace import Workspace
from nencarta.core.model_config import ModelConfig
from nencarta.core.enumerations import FloodMapMode


def prepare_inputs_for_dem(workspace: Workspace) -> ModelConfig:
    DEM = workspace.assigned_dem
    if not DEM.suffix.lower() in {'.tif', '.vrt', '.img', '.tiff', '.btiff', '.bil', '.dem'}:
        return ModelConfig(None, [], [], workspace.configs.mapper, workspace.configs.quiet)
    
    assign_and_validate_dem(workspace)
    Create_BaseLine_Manning_n_File_ESA(workspace)
    make_land_cover(workspace)
    make_stream_geometry(workspace)
    make_reanalysis_file(workspace)
    make_flood_flow_file_from_base_max_file(workspace, ['COMID', 'rp2'])

    if not workspace.configs.mapper.is_curve2flood_fldpln_mapper():
        make_stream_raster(workspace)

    if workspace.configs.clean_dem or not workspace.configs.disable_bathymetry:
        make_water_mask(workspace)

    if workspace.configs.clean_dem:
        define_configs_for_dem_cleaning(workspace)
        make_clean_dem(workspace)

    if workspace.configs.mapper.is_curve2flood_fldpln_mapper():
        make_fldpln_inputs(workspace)

    if workspace.configs.floodmap_mode == FloodMapMode.FORECAST:
        flow_files = make_flow_file_from_forecast(workspace)
        if workspace.configs.remove_old_forecast_files:
            remove_old_forecast_files(workspace, flow_files[0].name)
    elif workspace.configs.floodmap_mode == FloodMapMode.USER:
        flow_files = assign_user_flow_files(workspace.configs)
    elif workspace.configs.floodmap_mode == FloodMapMode.RETURN_PERIOD:
        flow_files = make_return_period_flow_file(workspace)
    else:
        flow_files = []

    arc_config = define_arc_configs(workspace)

    mapper_configs = []
    fist_inputs = []
    for flow_file in flow_files:
        mapper_configs.append(define_mapper_configs(workspace, flow_file))
        if workspace.configs.make_fist_inputs:
            fist_inputs.extend(get_fist_inputs(workspace, flow_file))

    return ModelConfig(
        arc_config,
        mapper_configs, 
        fist_inputs,
        workspace.configs.mapper,
        workspace.configs.quiet)