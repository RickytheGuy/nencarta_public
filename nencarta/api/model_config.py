from pathlib import Path
from functools import cached_property

from nencarta.api.inspect_config import inspect_config
from nencarta.api.enumerations import Mapper

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

    @cached_property
    def vdt_exists(self) -> bool:
        return self.arc_config.exists() and Path(inspect_config(self.arc_config, "Print_VDT_Database")).exists()