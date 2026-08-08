# build-in imports
import os
import json
import pprint
import argparse
from pathlib import Path
import multiprocessing as mp
from contextlib import redirect_stdout, redirect_stderr

# local imports
from nencarta.logger import LOG
from nencarta.workspace import Workspace
from nencarta.pipeline import build_pipeline
from nencarta.core.enumerations import Mapper
from nencarta.core.configs import NencartaConfig
from nencarta.Download_Process_ForecastData import Download_USGS_DEM_Data_Using_WarningFlag_Data


def _resolve_parallel_settings(data, cli_parallel=None, cli_num_workers=None):
    """
    Resolve (parallel: bool, num_workers: int) from JSON content and optional
    CLI overrides. CLI overrides win when provided. Defaults to serial.
    """
    # Support both 'parallel' and legacy 'run_parallel' in JSON
    json_parallel = data.get("parallel", data.get("run_parallel"))
    json_workers = data.get("num_workers", data.get("workers"))

    # Choose final parallel flag
    parallel = cli_parallel if cli_parallel is not None else json_parallel
    if parallel is None:
        parallel = False  # safer, deterministic default

    # Choose final worker count
    num_workers = cli_num_workers if cli_num_workers is not None else json_workers
    if parallel:
        try:
            num_workers = int(num_workers) if num_workers is not None else (os.cpu_count() or 1)
        except (TypeError, ValueError):
            num_workers = os.cpu_count() or 1
        num_workers = max(1, num_workers)
    else:
        num_workers = 1

    return bool(parallel), num_workers
    
def process_json_input_serial(json_file):
    """Process input from a JSON file."""
    with open(json_file, 'r') as file:
        LOG.info(f'Opening JSON file: {json_file}')
        data: dict = json.load(file)
        LOG.info(f"{os.linesep}{pprint.pformat(data)}")
    
    watersheds: list[dict] = data.get("watersheds", [])
    timer = Timer()
    for watershed in watersheds:
        process_watershed(watershed, timer)

    return timer

def process_single_watershed(watershed: dict):
    """Run a single watershed processing job and log to a file."""
    watershed_name = watershed.get("name", "unnamed")
    safe_name = watershed_name.replace(" ", "_").replace("/", "_")
    log_dir = os.path.join("logs")
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, f"{safe_name}.log")

    with open(log_path, 'w') as log_file, redirect_stdout(log_file), redirect_stderr(log_file):
        try:
            timer = Timer()
            process_watershed(watershed, timer)
            return timer
        except Exception as e:
            LOG.error(f"Exception during processing {watershed_name}: {e}")

def _resolve_workspace_dems(configs: NencartaConfig) -> list:
    """Return DEM inputs that should become workspaces."""
    if configs.dem_dir:
        if configs.use_warning_flags_to_download_dem:
            if configs.forensic_forecast_date:
                return Download_USGS_DEM_Data_Using_WarningFlag_Data(
                    configs.geoglows_vpu,
                    configs.dem_dir,
                    configs.forensic_forecast_date,
                )
            return []
        return list(Path(configs.dem_dir).glob(configs.dem_filter))

    if configs.dem:
        return [configs.dem]

    return []

def process_watershed(input_dict: dict):
    """Process one watershed dictionary through the workspace pipeline.

    The ``timer`` argument is kept for compatibility with older callers. The
    current pipeline manages timing through task-level profiling.
    """
    configs = NencartaConfig(input_dict)
    dem_list = _resolve_workspace_dems(configs)
    
    if not dem_list:
        if configs.source_dems and configs.bbox:
            dem_list = [None]
        else:
            LOG.warning("No DEMs found in the specified folder and no source DEMs provided; cannot run pipeline.")
            return
        
    workspaces = [Workspace(configs, dem) for dem in dem_list]
    run_pipeline(workspaces)
    return

def run_pipeline(workspaces: list[Workspace], executor=None):
    profile = workspaces[0].configs.profile and not workspaces[0].configs.parallel
    parallel = workspaces[0].configs.parallel
    pipeline = build_pipeline(profile)

    # pipeline.visualize_graphviz(filename='test.png', hide_default_args=True, orient='TB')

    try:
        pipeline.map(
            {'workspace': workspaces},
            parallel=parallel,
            scheduling_strategy='eager',
            show_progress="rich",
            executor=executor,
            # error_handling='continue'
            )
        
        if profile:
            pipeline.print_profiling_stats()
    except Exception as e:
        if not pipeline.error_snapshot:
            raise e
    finally:
        if pipeline.error_snapshot:
            print(pipeline.error_snapshot.traceback)
            print(pipeline.error_snapshot)

    LOG.info(f"Finished processing")
    return

def process_json_input(json_file, parallel=None, num_workers=None):
    """
    Process watersheds defined in a JSON file, optionally in parallel.

    If 'parallel' / 'num_workers' are None, they are resolved from the JSON file
    via keys "parallel" (or legacy "run_parallel") and "num_workers" (or "workers").
    CLI overrides take precedence when provided.
    """
    with open(json_file, 'r') as file:
        LOG.info(f"Reading input file: {json_file}")
        data = json.load(file)

    watersheds = data.get("watersheds", [])
    if not watersheds:
        LOG.warning("No watersheds found.")
        return

    # Resolve final parallel settings (JSON + optional CLI override)
    parallel, num_workers = _resolve_parallel_settings(
        data,
        cli_parallel=parallel,
        cli_num_workers=num_workers
    )
    LOG.info(f"Parallel: {parallel} | Workers: {num_workers}")
    if parallel and len(watersheds) > 1:
        # Parallel path: one process per watershed, logging handled by process_single_watershed
        with mp.Pool(processes=num_workers) as pool:
            pool.imap_unordered(process_single_watershed, watersheds)
    else:
        # Serial path: reuse existing serial routine (retains your validation flow)
        process_json_input_serial(json_file)

def rename_cli_keys(input_dict: dict) -> dict:
    """Rename CLI argument keys to match watershed dictionary keys."""
    key_mapping = {
        "watershed":"name",
    }
    return {key_mapping.get(key, key): value for key, value in input_dict.items()}

def process_cli_arguments(args):
    """Process input from CLI arguments."""
    input_dict = vars(args)
    input_dict = rename_cli_keys(input_dict)
    process_watershed(input_dict)

class RequiredIfFloodWaterAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values)

        if namespace.flood_waterlc_and_strm_cells and values is None:
            parser.error("--land_watervalue is required when --flood_waterlc_and_strm_cells is set to True.")


def main():
    parser = argparse.ArgumentParser(description="Flood Mapping Script")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand for JSON input
    json_parser = subparsers.add_parser("json", help="Process watersheds from a JSON file")
    json_parser.add_argument("json_file", type=str, help="Path to the JSON file")
    # Let CLI override JSON for parallel choice
    grp = json_parser.add_mutually_exclusive_group()
    grp.add_argument("--parallel", dest="parallel", action="store_const", const=True, default=None,
                     help="Force parallel processing (overrides JSON).")
    grp.add_argument("--serial", dest="parallel", action="store_const", const=False,
                     help="Force serial processing (overrides JSON).")
    json_parser.add_argument("--num_workers", type=int, default=None,
                             help="Number of workers when running in parallel (overrides JSON).")

    # # Subcommand for JSON input
    # json_parser = subparsers.add_parser("json", help="Process watersheds from a JSON file")
    # json_parser.add_argument("json_file", type=str, help="Path to the JSON file")

    # Subcommand for CLI input
    cli_parser = subparsers.add_parser("cli", help="Process watershed parameters via CLI")
    cli_parser.add_argument("watershed", type=str, help="Watershed name")
    cli_parser.add_argument("flowline", type=str, help="Path to the flowline shapefile")
    cli_parser.add_argument("dem_dir", type=str, help="Directory containing DEM files")
    cli_parser.add_argument("output_dir", type=str, help="Directory where results will be saved")
    cli_parser.add_argument("--disable_bathymetry", action="store_true", help="Disable ARC, Curve2Flood, and FloodSpreader bathymetry estimation inputs")
    cli_parser.add_argument("--bathy_use_banks", action="store_true", help="Use bathy banks for processing")
    cli_parser.add_argument("--flood_waterlc_and_strm_cells", action="store_true",
                        help="In the flood inundation maps it shows water related land use and stream raster cells as flooded")
    cli_parser.add_argument("--land_watervalue", type=int, action=RequiredIfFloodWaterAction,
                        help="Land water value in the land cover raster (Required if --flood_waterlc_and_strm_cells is True)")
    cli_parser.add_argument("--clean_dem", action="store_true", help="Clean DEM data before processing")
    cli_parser.add_argument("--process_stream_network", action="store_true", help="Clean DEM data before processing")
    cli_parser.add_argument("--mapper", type=str, default="Curve2Flood-Kernel Weighted", choices=Mapper.list_names(), help="Mapping method")
    cli_parser.add_argument("--use_specified_depth_for_bathy_mask", action="store_true", help="Specify a depth for FloodSprederPy to use for bathymetry masking")
    cli_parser.add_argument("--age_of_forecast_days", type=int, default=7, help="Age of forecast in days")
    cli_parser.add_argument("--find_banks_based_on_landcover", action="store_true", help="Use landcover data for finding banks when estimating bathymetry")
    cli_parser.add_argument("--specify_depths_for_bathy_mask", type=float, nargs=2, help="Specify two floats for bathymetry depth mask when '--use_specified_depth_for_bathy_mask' is True")
    cli_parser.add_argument("--create_reach_average_curve_file", action="store_true", help="Create a reach average curve file instead of one that varies for each stream cell")
    cli_parser.add_argument("--forensic_forecast_date", type=str, default=None, help="Forensic forecast date in YYYYMMDD format (defaults to most recent forecast) unless this argument is provided")
    cli_parser.add_argument("--forensic_forecast_hour", type=int, default=0, choices=[i for i in range (0,24)], help="Forensic forecast hour (defaults to 0) unless this argument is provided")
    cli_parser.add_argument("--specified_bathyflow_field", type=str, default="p_exceed_50", help="Specify the streamflow field in the streamflow reanalysis file that will be used for bathymetry estimation  (defaults to 'p_exceed_50') ")
    cli_parser.add_argument("--specified_highflow_field", type=str, default="rp100_premium", help="Specify the highflow field in the streamflow reanalysis file that will be used by ARC as the highest flow for the VDT database and curvefile creation (defaults to 'rp100_premium')")
    cli_parser.add_argument("--StrmOrder_Field", type=str, default=None, help="Stream order field in the stream shapefile (optional)")
    cli_parser.add_argument("--StrmOrder_Lower", type=int, default=None, help="Upper bound for stream order (optional)")
    cli_parser.add_argument("--StrmOrder_Upper", type=int, default=None, help="Upper bound for stream order (optional)")
    cli_parser.add_argument("--q_baseflow_threshold", type=float, default=None, help="Drop streams whose baseflow is below this threshold (optional)")
    cli_parser.add_argument("--use_warning_flags_to_download_dem", action="store_true", help="Use warning flags to download DEM data")
    cli_parser.add_argument("--geoglows_vpu", type=int, default=None, help="GEOGLOWS VPU ID (required if --use_warning_flags_to_download_dem is set to True)")
    cli_parser.add_argument("--lake_filter_json", type=str, default=None, help="Path to the lake filter JSON file (optional)")
    cli_parser.add_argument("--estimate_consequences", action="store_true", help="Estimate consequences using go-consequences")
    cli_parser.add_argument("--streamflow_source", type=str, default="GEOGLOWS", choices=["NWM", "GEOGLOWS"], help="Streamflow source for NenCarta (defaults to GEOGLOWS)")
    cli_parser.add_argument("--nwm_api_key", type=str, default=None, help="NWM API key (required when --streamflow_source is NWM)")
    cli_parser.add_argument("--overwrite_floodmaps", action="store_true", help="Overwrite existing forecast flood maps")
    cli_parser.add_argument("--remove_old_forecast_files", action="store_true", help="Remove old forecast files before processing")
    cli_parser.add_argument("--make_fist_inputs", action="store_true", help="Make FIST inputs after processing")
    cli_parser.add_argument("--move_stream_network_to_thalweg", action="store_true", help="Move stream network to the thalweg of the DEM")
    cli_parser.add_argument("--new_strm_threshold_km2", type=float, default=None, help="The stream threshold for creating a new stream network for the DEM that you will be using. Use in conjunction with move_stream_network_to_thalweg and Curve2Flood-FLDPLNpy")
    cli_parser.add_argument("--min_match_score", type=float, default=None, help="The score needed to conflate the new DEM based network with the one provided by the as the `flowline` input")
    gui_parser = subparsers.add_parser("gui", help="Summon the GUI application")

    args = parser.parse_args()

    

    if args.command == "json":
        LOG.info('Processing ' + str(args.json_file))
        process_json_input(args.json_file, parallel=args.parallel, num_workers=args.num_workers)
    elif args.command == "cli":
        process_cli_arguments(args)
    elif args.command == "gui":
        from nencarta.gui_app import run_gui
        run_gui()


if __name__ == "__main__":
    main()
