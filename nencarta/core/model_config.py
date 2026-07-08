from pathlib import Path

from nencarta.core.inspect_config import inspect_config
from nencarta.core.enumerations import Mapper

class ModelConfig():
    def __init__(
            self, 
            arc_config: Path, 
            mapper_configs: list[Path],
            fist_inputs: list[tuple],
            mapper: Mapper,
            quiet: bool = False):
        self.arc_config = arc_config
        self.mapper_configs = mapper_configs
        self.fist_inputs = fist_inputs
        self.mapper = mapper
        self.quiet = quiet

    @property
    def vdt_exists(self) -> bool:
        return self.arc_config.exists() and Path(inspect_config(self.arc_config, "Print_VDT_Database")).exists()
    
    def __repr__(self) -> str:
        return f"ModelConfig(arc_config={self.arc_config}, mapper={self.mapper})"