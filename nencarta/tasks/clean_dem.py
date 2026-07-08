from pathlib import Path

from nencarta.logger import LOG
from nencarta import DEM_Cleaner
from nencarta.workspace import Workspace
from nencarta.core.configs import NencartaConfig
from nencarta.tasks.run_models import _run_mapper, run_arc_for_initial_floodmap
from nencarta.tasks.make_flow_file_from_reanalysis import _make_flood_flow_file_from_base_max_file

def make_clean_dem(workspace: Workspace) -> Path:
    # Run the DEM Cleaner Program, if you wanna
    configs = workspace.configs
    if workspace.DEM_File_Clean.exists() and not configs.process_stream_network:
        return workspace.DEM_File_Clean

    if configs.process_stream_network:
        # Ensure we regenerate ARC + DEM cleaner products when overwriting, so the
        # outputs match a fresh legacy run (important for ignore/*_original fixtures).
        for path in (
            workspace.DEM_File_Clean,
            workspace.Curve_File_Initial,
            workspace.VDT_File_Initial,
        ):
            path.unlink(missing_ok=True)

        dem_name = workspace.DEM_File_Clean.name
        for intermediate_name in (
            f"Elev_Streams_{dem_name}",
            f"FLOOD_IMPACT_{dem_name}",
            f"SEED_CONNECT_{dem_name}",
        ):
            intermediate_path = workspace.dem_updated_folder / intermediate_name
            intermediate_path.unlink(missing_ok=True)
        
    if not workspace.Curve_File_Initial.exists():
        run_arc_for_initial_floodmap(workspace.ARC_FileName_for_DEM_Cleaner, configs)
        if not workspace.Curve_File_Initial.exists():
            raise FileNotFoundError(f"Expected Curve2Flood output file not found at {workspace.Curve_File_Initial} after running ARC for initial floodmap.")
    else:
        LOG.info(f"{workspace.Curve_File_Initial} exists and we aren't making it again...")

    if configs.use_specified_depth_for_bathy_mask:
        LOG.info(f"Executing Curve2Flood using {workspace.ARC_FileName_for_DEM_Cleaner}")
        _run_mapper(workspace.ARC_FileName_for_DEM_Cleaner, configs)
    
    OutputID = 'COMID'
    Q_Fraction = 0.10
    TopWidthPlausibleLimit = 600
    search_dist_for_min_elev = 10
    search_dist_perp_cells = 10 # this was 40
    FlowFileName = workspace.FLOW_Folder / f"{workspace.FileName}_Flow_COMID_Q.txt"
    _make_flood_flow_file_from_base_max_file(
        reanalysis_flow_file=workspace.DEM_Reanalsyis_FlowFile,
        out_file=FlowFileName,
        columns=[OutputID, 'p_exceed_50'],
    )

    workspace.dem_updated_folder.mkdir(parents=True, exist_ok=True)

    # start time for the simulation
    DEM_Cleaner.DEM_Cleaner_Program(OutputID, 
                                    workspace.DEM_StrmShp, 
                                    workspace.assigned_dem.parent, 
                                    [workspace.assigned_dem.name], 
                                    [workspace.STRM_File_Clean], 
                                    workspace.dem_updated_folder, 
                                    FlowFileName, 
                                    workspace.Curve_File_Initial, 
                                    workspace.bathy_water_mask, 
                                    Q_Fraction, 
                                    TopWidthPlausibleLimit, 
                                    search_dist_for_min_elev, 
                                    search_dist_perp_cells)
    
    workspace.assigned_dem = workspace.DEM_File_Clean