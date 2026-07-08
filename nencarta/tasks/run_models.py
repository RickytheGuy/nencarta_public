import subprocess
from pathlib import Path

from nencarta.logger import LOG
from nencarta.workspace import Workspace
from nencarta.core.enumerations import Mapper
from nencarta.core.configs import NencartaConfig
from nencarta._constants import FLOODSPREADER_PATH
from nencarta.core.floodmapper_output import FloodMapperBulkOutput, FloodMapperOutput
from nencarta.core.model_config import ModelConfig

def _run_arc(config: Path, model_config: ModelConfig) -> None:
    from arc import Arc # Lazy import helps load faster
    Arc(str(config), quiet=model_config.quiet).run()

def _run_mapper(config_file: Path, model_config: ModelConfig) -> None:
    if model_config.mapper == Mapper.FLOODSPREADER:
        if not FLOODSPREADER_PATH.exists():
            raise FileNotFoundError(f"FloodSpreader script not found at {FLOODSPREADER_PATH}. Please ensure it is included in the nencarta package.")
        call_mapper = f'python "{FLOODSPREADER_PATH}" {config_file}'
        subprocess.call(call_mapper, shell=True)
    else:
        from curve2flood import Curve2Flood_MainFunction
        Curve2Flood_MainFunction(str(config_file), quiet=model_config.quiet)

def run_arc_bathymetry(model_config: ModelConfig) -> ModelConfig:
    LOG.info("Running ARC...")
    _run_arc(model_config.arc_config, model_config)
    
    return model_config

def run_mapper_bathymetry(model_config: ModelConfig, workspace: Workspace) -> ModelConfig:
    if not model_config.vdt_exists:
        return model_config
    
    if not workspace.configs.disable_bathymetry:
        LOG.info("Running flood mapper to generate bathymetry...")
        _run_mapper(model_config.arc_config, model_config)

    return model_config

def run_mapper_floodmaps(model_config: ModelConfig) -> FloodMapperBulkOutput:
    LOG.info("Running flood mapper to generate flood maps...")
    for config in model_config.mapper_configs:
        if model_config.vdt_exists:
            _run_mapper(config, model_config)

    return FloodMapperBulkOutput(model_config.mapper_configs)

def run_arc_for_initial_floodmap(config: Path,
            configs: NencartaConfig) -> Path:
    LOG.info("Running ARC to generate initial floodmap for DEM cleaner...")
    _run_arc(config, configs)
    
    return config