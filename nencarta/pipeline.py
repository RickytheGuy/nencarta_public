from pipefunc import Pipeline, PipeFunc

import nencarta.tasks as tasks

# def get_pipeline_names(configs: NencartaConfig) -> tuple[str, str, bool]:
#     move_streams = (
#         configs[0].move_stream_network_to_new_locations
#         or configs[0].mapper.is_curve2flood_fldpln_mapper()
#     )

#     if configs[0].clean_dem:
#         return "intermediate_dem", "final_stream_geometry", move_streams

#     if move_streams:
#         return "final_dem", "intermediate_stream_geometry", move_streams

#     return "final_dem", "final_stream_geometry", move_streams

# def add_base_inputs(funcs: list, dem_name, stream_name):
#     funcs.extend([
#         PipeFunc(
#             tasks.Create_BaseLine_Manning_n_File_ESA,
#             output_name="mannings_n_file",
#             mapspec="workspace[i] -> mannings_n_file[i]",      
#             ),
#         PipeFunc(
#             tasks.assign_and_validate_dem, 
#             dem_name,
#             mapspec=f"workspace[i], configs[i] -> {dem_name}[i]",
#         ),
#         PipeFunc(
#             tasks.make_land_cover,
#             "land_cover",
#             renames={'dem': dem_name},
#             mapspec=f"workspace[i], {dem_name}[i], configs[i] -> land_cover[i]",
#         ),
#         PipeFunc(
#             tasks.make_stream_geometry,
#             stream_name,
#             renames={'dem': dem_name},
#             mapspec=f"workspace[i], {dem_name}[i], configs[i] -> {stream_name}[i]",
#         ),
#     ])

# def add_stream_relocation(funcs: list, dem_name: str):
#     funcs.append(
#         PipeFunc(
#             tasks.make_fldpln_inputs,
#             'final_stream_geometry',
#             renames={'dem': dem_name, 'stream_geometry': 'intermediate_stream_geometry'},
#             mapspec=f"workspace[i], {dem_name}[i], intermediate_stream_geometry[i], configs[i] -> final_stream_geometry[i]",
#         )
#     )   

# def add_clean_dem(funcs: list, dem_name: str):
#     funcs.extend([
#         PipeFunc(
#             tasks.define_configs_for_dem_cleaning,
#             'clean_dem_config',
#             renames={'dem': dem_name, 'stream_geometry': 'final_stream_geometry'},
#             mapspec=f"workspace[i], {dem_name}[i], final_stream_geometry[i],  mannings_n_file[i], reanalysis_flow_file[i], stream_raster[i], two_yr_flow_file[i], land_cover[i], configs[i] -> clean_dem_config[i]",
#         ),
#         PipeFunc(
#             tasks.make_clean_dem,
#             "final_dem",
#             renames={'dem': dem_name, 'stream_geometry': 'final_stream_geometry'},
#             mapspec=f"workspace[i], {dem_name}[i], final_stream_geometry[i], clean_dem_config[i], reanalysis_flow_file[i], bathy_water_mask[i], stream_raster[i], configs[i] -> final_dem[i]",
#         ),
#     ])


# def add_bathymetry(funcs: list, configs: NencartaConfig, dem_name: str):
#     ...

# def add_flow_preparation(funcs: list, configs: NencartaConfig, dem_name: str):
#     funcs.extend([
#         PipeFunc(
#             tasks.make_reanalysis_file,
#             "reanalysis_flow_file",
#             renames={'stream_geometry': 'final_stream_geometry'},
#             mapspec=f"workspace[i], final_stream_geometry[i], configs[i] -> reanalysis_flow_file[i]",
#         ),
#         PipeFunc(
#             tasks.make_stream_raster,
#             "stream_raster",
#             renames={'dem': dem_name, 'stream_geometry': 'final_stream_geometry'},
#             mapspec=f"workspace[i], {dem_name}[i], final_stream_geometry[i], configs[i] -> stream_raster[i]",
#         ),
#         PipeFunc(
#             tasks.make_flood_flow_file_from_base_max_file,
#             "two_yr_flow_file",
#             defaults={'columns': ['COMID', 'rp2'],},
#             mapspec=f"workspace[i], reanalysis_flow_file[i], configs[i] -> two_yr_flow_file[i]",
#         )
#     ])



# def add_flow_generation(funcs: list, configs: NencartaConfig):
#     if configs[0].floodmap_mode == FloodMapMode.FORECAST:
#         funcs.append(
#             PipeFunc(
#                 tasks.make_flow_file_from_forecast,
#                 "flow_files",
#                 renames={'stream_geometry': 'final_stream_geometry'},
#                 mapspec='workspace[i], final_stream_geometry[i], configs[i] -> flow_files[i, j]',
#             )
#         )
#     elif configs[0].floodmap_mode == FloodMapMode.USER:
#         funcs.append(
#             PipeFunc(
#                 tasks.assign_user_flow_files,
#                 "flow_files",
#                 mapspec='configs[i] -> flow_files[i, j]',
#             )
#         )
#     elif configs[0].floodmap_mode == FloodMapMode.RETURN_PERIOD:
#         funcs.append(
#             PipeFunc(
#                 tasks.make_return_period_flow_file,
#                 "flow_files",
#                 mapspec='workspace[i], reanalysis_flow_file[i], configs[i] -> flow_files[i, j]',
#             )
#         )
#     else:
#         LOG.error(f"Invalid floodmap_mode: {configs[0].floodmap_mode}")
#         raise ValueError(f"Invalid floodmap_mode: {configs[0].floodmap_mode}")

# def add_mapper_configs(funcs: list):
#     funcs.extend([
#         PipeFunc(
#             tasks.define_arc_configs,
#             'arc_bathy_config',
#             renames={'dem': 'final_dem', 'stream_geometry': 'final_stream_geometry'},
#             mapspec="workspace[i], mannings_n_file[i], stream_raster[i], final_stream_geometry[i], final_dem[i], reanalysis_flow_file[i], two_yr_flow_file[i], bathy_water_mask[i], land_cover[i], configs[i] -> arc_bathy_config[i]",
#         ),
#         PipeFunc(
#             tasks.define_mapper_configs,
#             "c2f_configs",
#             renames={'dem': 'final_dem', 'stream_geometry': 'final_stream_geometry', 'flow_file': 'flow_files'},
#             mapspec="workspace[i], mannings_n_file[i], flow_files[i, j], final_stream_geometry[i], stream_raster[i], final_dem[i], land_cover[i], configs[i] -> c2f_configs[i, j]"
#         ),
#         PipeFunc(
#             tasks.run_arc_bathymetry,
#             "mapper_bathy_config",
#             renames={'config': 'arc_bathy_config'},
#             mapspec="arc_bathy_config[i], configs[i] -> mapper_bathy_config[i]",
#         ),
#         PipeFunc(
#             tasks.run_mapper_bathymetry,
#             "vdt",
#             renames={'config': 'mapper_bathy_config'},
#             mapspec="mapper_bathy_config[i], workspace[i], configs[i] -> vdt[i]",
#         ),
#         PipeFunc(
#             tasks.run_mapper_floodmaps,
#             "mapper_outputs",
#             renames={'config': 'c2f_configs'},
#             mapspec="c2f_configs[i, j], configs[i] -> mapper_outputs[i, j]",
#         ),
#         PipeFunc(
#             tasks.unbuffer_maps,
#             "unbuffered_floodmaps",
#             renames={'floodmapper_output': 'mapper_outputs'},
#             mapspec="mapper_outputs[i, j], workspace[i], configs[i] -> unbuffered_floodmaps[i, j]",
#         )
#     ])

# def add_optional_outputs(funcs: list, configs: NencartaConfig):
#     if configs[0].clean_dem or not configs[0].disable_bathymetry:
#         funcs.append(
#             PipeFunc(
#                 tasks.make_water_mask,
#                 "bathy_water_mask",
#                 mapspec=f"workspace[i], stream_raster[i], land_cover[i], configs[i] -> bathy_water_mask[i]",
#             )
#         )

#     if configs[0].floodmap_mode == FloodMapMode.FORECAST and configs[0].remove_old_forecast_files:
#         funcs.append(
#             PipeFunc(
#                 tasks.remove_old_forecast_files,
#                 'remove_old_forecast_files_done',
#                 mapspec='workspace[i], configs[i] -> remove_old_forecast_files_done[i]',
#             )
#         )

#     if configs[0].make_fist_inputs:
#         funcs.extend([
#             PipeFunc(
#                 tasks.get_fist_inputs,
#                 "fist_inputs",
#                 renames={'stream_geometry': 'final_stream_geometry', 'flow_file': 'flow_files'},
#                 mapspec="flow_files[j], configs[i] -> fist_inputs[j]",
#             ),
#             PipeFunc(
#                 tasks.run_fist,
#                 "fist_outputs",
#                 renames={'args': 'fist_inputs'},
#                 mapspec="fist_inputs[j] -> fist_outputs[j]",
#             )
#         ])

#     if configs[0].estimate_consequences:
#         funcs.extend([
#             PipeFunc(
#                 tasks.get_consequences_tasks,
#                 "consequences_tasks",
#                 renames={'floodmapper_output': 'mapper_outputs'},
#                 mapspec="mapper_outputs[i] -> consequences_tasks[i]",
#             ),
#             PipeFunc(
#                 tasks.run_consequences,
#                 "run_consequences",
#                 renames={'docker_command': 'consequences_tasks'},
#                 mapspec="consequences_tasks[i] -> run_consequences[i]",
#             )
#         ])

def build_pipeline(profile: bool) -> Pipeline:
    funcs = [
        PipeFunc(
            tasks.prepare_inputs_for_dem,
            "model_configs",
            mapspec="workspace[i] -> model_configs[i]",
        ),
        PipeFunc(
            tasks.run_arc_bathymetry,
            "mapper_arc_config",
            renames={'model_config': 'model_configs'},
            mapspec="model_configs[i], workspace[i] -> mapper_arc_config[i]",
        ),
        PipeFunc(
            tasks.run_mapper_bathymetry,
            "mapper_bathy_config",
            renames={'model_config': 'mapper_arc_config'},
            mapspec="mapper_arc_config[i], workspace[i] -> mapper_bathy_config[i]",
        ),
        PipeFunc(
            tasks.run_fldpln_library,
            "mapper_ready_config",
            renames={'model_config': 'mapper_bathy_config'},
            mapspec="mapper_bathy_config[i], workspace[i] -> mapper_ready_config[i]",
        ),
        PipeFunc(
            tasks.run_mapper_floodmaps,
            "mapper_outputs",
            renames={'model_config': 'mapper_ready_config'},
            mapspec="mapper_ready_config[i], workspace[i] -> mapper_outputs[i]",
        ),
        PipeFunc(
            tasks.unbuffer_maps,
            "unbuffered_floodmaps",
            renames={'floodmapper_output': 'mapper_outputs'},
            mapspec="mapper_outputs[i], workspace[i] -> unbuffered_floodmaps[i]",
        ),
        PipeFunc(
            tasks.run_fist,
            "fist_outputs",
            renames={'model_config': 'model_configs', 'flag': 'mapper_arc_config'}, # Flag to indicate we need to run ARC first before FIST
            mapspec="model_configs[i], mapper_arc_config[i] -> fist_outputs[i]",
        ),
        PipeFunc(
            tasks.run_consequences,
            "consequences_outputs",
            renames={'floodmapper_bulk_output': 'mapper_outputs'},
            mapspec="mapper_outputs[i], workspace[i] -> consequences_outputs[i]",
        )
    ]

    return Pipeline(
        funcs,
        profile=profile
    )