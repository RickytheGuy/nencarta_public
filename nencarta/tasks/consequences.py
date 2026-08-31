import json
import subprocess

from osgeo import gdal

from nencarta.logger import LOG
from nencarta.core.floodmapper_output import FloodMapperBulkOutput
from nencarta.workspace import Workspace

gdal.UseExceptions()

def _Create_Go_Consequence_GeoJSON(Consequences_JSON_Path, Container_Hazard_Path, Consequences_Output_GPKG_File):
    config = {
        "structure_provider_info": {
            "structure_provider_type": "NSIAPI"
        },
        "hazard_provider_info": {
            "hazards": [
                {
                    "hazard_parameter_type": "depth",
                    "hazard_provider_file_path": Container_Hazard_Path
                }
            ]
        },
        "results_writer_info": {
            "results_writer_type": "GPKG",
            "output_file_path": f"/data/Consequences/{Consequences_Output_GPKG_File}"
        }
    }

    with open(Consequences_JSON_Path, "w", encoding="utf-8") as f:
        json.dump(config, f, indent=4)

    return

def _get_consequences_tasks(floodmapper_bulk_output: FloodMapperBulkOutput, workspace: Workspace) -> list:
    if not workspace.configs.estimate_consequences:
        return []

    workspace.Consequences_Folder.mkdir(parents=True, exist_ok=True)
    
    tasks = []
    for floodmapper_output in floodmapper_bulk_output.outputs:
        depth_file = floodmapper_output.depth_map
        if depth_file is None:
            LOG.warning("No depth file found in floodmapper output, skipping consequences estimation")
            continue
        
        # Convert the depth raster to WGS84 coordinate system before processing
        depth_file_wgs84 = depth_file.with_name(depth_file.stem + '_WGS84' + depth_file.suffix)
        LOG.info(f"Converting depth raster to WGS84: {depth_file} -> {depth_file_wgs84}")
        try:
            warp_options = gdal.WarpOptions(
                format="GTiff",
                dstSRS="EPSG:4326",
                resampleAlg=gdal.GRA_Bilinear,
                creationOptions=["COMPRESS=DEFLATE"],
            )
            ds = gdal.Warp(depth_file_wgs84, depth_file, options=warp_options)
            if ds is None:
                LOG.error(f"Failed to convert {depth_file} to WGS84, using original file")
                depth_file_wgs84 = depth_file
            else:
                ds = None
                LOG.info(f"Successfully converted to WGS84: {depth_file_wgs84}")
        except Exception as e:
            LOG.error(f"Error converting to WGS84: {e}, using original file")
            depth_file_wgs84 = depth_file
        
        # Resolve container path based on workspace output mount (/data)
        rel_depth_path = depth_file_wgs84.relative_to(workspace.output_dir).as_posix()
        Container_Hazard_Path = f"/data/{rel_depth_path}"

        Forecast_Flood_Depth_Raster_Name = depth_file_wgs84.name
        Consequences_JSON_File = Forecast_Flood_Depth_Raster_Name.replace('.tif', '_consequences.json') 
        Consequences_JSON_Path = workspace.Consequences_Folder / Consequences_JSON_File
        Consequences_Output_GPKG_File = Consequences_JSON_File.replace('.json', '.gpkg')

        LOG.info(f"Creating consequences file {Consequences_JSON_Path}")
        _Create_Go_Consequence_GeoJSON(Consequences_JSON_Path, Container_Hazard_Path, Consequences_Output_GPKG_File)
        # run the go-consequences Docker container
        source_dir = workspace.output_dir
        source_dir = source_dir.resolve().as_posix()

        docker_command = [
            "docker",
            "run",
            "--rm",
            "--mount",
            f"type=bind,source={source_dir},target=/data",
            "go-consequences:latest",
            f"/data/Consequences/{Consequences_JSON_File}",
        ]
        tasks.append(docker_command)

    return tasks

def run_consequences(floodmapper_bulk_output: FloodMapperBulkOutput, workspace: Workspace) -> None:
    consequences_tasks = _get_consequences_tasks(floodmapper_bulk_output, workspace)
    for docker_command in consequences_tasks:
        LOG.info(f"Running consequences estimation with command: {' '.join(docker_command)}")
        try:
            result = subprocess.run(docker_command, check=True, capture_output=True, text=True)
            LOG.info(f"Consequences estimation completed successfully: {result.stdout}")
        except subprocess.CalledProcessError as e:
            LOG.error(f"Error running consequences estimation: {e.stderr}")