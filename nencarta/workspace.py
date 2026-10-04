from pathlib import Path

from nencarta.core.enumerations import  Mapper, StreamflowSource
from nencarta.core.configs import NencartaConfig

# whitebox_tools splits its arguments on '=' and strips quotes out of them, so it cannot open a
# path containing any of these: '=' makes it panic, and quotes make it report a missing file
WHITEBOX_UNSAFE_PATH_CHARACTERS = ('=', "'", '"')

# The folders a workspace writes to, which 'folder_paths' can move, and where they go by default
WORKSPACE_FOLDERS = ('DEM', 'ARC_InputFiles', 'FloodMap', 'Bathymetry', 'DEM_Updated', 'STRM', 'LAND',
                     'FLOW', 'VDT', 'ESA_LC', 'FIST', 'Consequences', 'FlowDirection', 'FLDPLN')
DEFAULT_FOLDER_PATH = '{name}/{folder}'

# What each kind of file is called after its '<source>_<DEM>_' or '<DEM>_' prefix, which 'file_names' can change
FILE_BLURBS = (
    'dem', 'fixed', 'filled', 'fixed_Clean', 'Clean',
    'flowdir', 'flowacc', 'wtbx_derived',
    'StrmShp', 'STRM_Raster', 'STRM_Raster_Clean', 'matched', 'lakes', 'stream_info', 'fldpln_library',
    'LAND_Raster', 'AR_Manning_n_MED',
    'Reanalysis', '2yr_flow_initial', 'rp', 'forecast', 'Flow_COMID_Q',
    'ARC_Input_Bathy', 'ARC_Input_InitialFlood', 'ARC_Input_FloodForecast', 'ARC_Input',
    'VDT_Database_Initial', 'VDT_Database_Bathy', 'AP_Database_Bathy', 'CurveFile', 'CurveFile_Initial',
    'CurveFile_Bathy', 'XS', 'Representative_XS',
    'water_mask', 'ARC_Bathy', 'FS_Bathy',
    'ARC_Flood', 'ARC_Flood_Initial', 'ARC_Flood_Bathy', 'ARC_Depth', 'ARC_FloodDepth', 'ARC_FloodWSE', 'ARC_FloodVEL',
    'Seed', 'FIST',
)

def check_workspaces_do_not_share_folders(workspaces: list['Workspace']) -> None:
    """
    Raise if two workspaces with short_file_names would write to the same folder: without the DEM in
    their file names, their files would overwrite each other.
    """
    owners = {}
    for workspace in workspaces:
        if not workspace.configs.short_file_names:
            continue
        for folder in WORKSPACE_FOLDERS:
            path = workspace.folder_path(folder)
            other = owners.setdefault(path, workspace)
            if other is not workspace:
                raise ValueError(
                    f"With 'short_file_names', {other.FileName} and {workspace.FileName} would both write to {path}. "
                    "Put {dem} in 'folder_paths', e.g. {\"default\": \"{name}/{dem}/{folder}\"}, so each DEM has its own folders."
                )

class Workspace:
    """Paths and per-DEM state for one watershed processing run."""

    def __init__(self, configs: NencartaConfig, dem: Path | None = None):
        self.configs = configs
        self.watershed: str = configs.watershed_name
        self.output_dir = Path(configs.output_dir) / self.watershed
        self.mapper: Mapper = configs.mapper

        unknown_folders = set(configs.folder_paths or {}) - set(WORKSPACE_FOLDERS) - {'default'}
        if unknown_folders:
            raise ValueError(
                f"Watershed '{self.watershed}': 'folder_paths' has unknown folder(s) {sorted(unknown_folders)}. "
                f"Use 'default' or any of {list(WORKSPACE_FOLDERS)}."
            )

        file_names = configs.file_names or {}
        unknown_blurbs = set(file_names) - set(FILE_BLURBS)
        if unknown_blurbs:
            raise ValueError(
                f"Watershed '{self.watershed}': 'file_names' has unknown file(s) {sorted(unknown_blurbs)}. "
                f"Use any of {list(FILE_BLURBS)}."
            )
        bad_names = {blurb: name for blurb, name in file_names.items() if not isinstance(name, str) or any(sep in name for sep in '/\\')}
        if bad_names:
            raise ValueError(f"Watershed '{self.watershed}': 'file_names' must give each file a name with no '/' or '\\': {bad_names}.")

        self.strm_source = 'NWM' if configs.streamflow_source.is_nwm() else 'GEOGLOWS'
        self.setup_dem_name(dem)
        self.DEM_folder = self.folder_path('DEM')
        self.ARC_Folder = self.folder_path('ARC_InputFiles')
        self.flood_folder = self.folder_path('FloodMap')
        self.bathy_file_folder = self.folder_path('Bathymetry')
        self.dem_updated_folder = self.folder_path('DEM_Updated')
        self.strm_folder = self.folder_path('STRM')
        self.land_folder = self.folder_path('LAND')
        self.FLOW_Folder = self.folder_path('FLOW')
        self.VDT_Folder = self.folder_path('VDT')
        self.ESA_LC_Folder = self.folder_path('ESA_LC')
        self.FIST_Folder = self.folder_path('FIST')
        self.Consequences_Folder = self.folder_path('Consequences')
        self.Flow_Direction_Folder = self.folder_path('FlowDirection')
        self.FLDPLN_Folder = self.folder_path('FLDPLN')
        
        configs_mannings_n = configs.mannings_text_file
        if configs_mannings_n:
            configs_mannings_n_path = Path(configs_mannings_n)
            if not configs_mannings_n_path.is_file():
                raise FileNotFoundError(f"Provided Manning's n text file not found: {configs_mannings_n}")
            self.mannings_n_text_file = configs_mannings_n_path
        else:
            self.mannings_n_text_file = self.file_path(self.land_folder, 'AR_Manning_n_MED', 'txt', None)

        self.floodmap_mode = configs.floodmap_mode

        self.setup_dem(dem)
        dem_updated, flow_direction = self.dem_updated_folder, self.Flow_Direction_Folder
        self.fixed_dem = self.file_path(dem_updated, 'fixed', 'tif')
        self.filled_dem = self.file_path(dem_updated, 'filled', 'tif')
        self.flowdir = self.file_path(flow_direction, 'flowdir', 'tif')
        self.flowacc = self.file_path(flow_direction, 'flowacc', 'tif')
        self.new_StrmShp = self.file_path(flow_direction, 'wtbx_derived', 'shp')
        self.whitebox_stream_raster = self.file_path(flow_direction, 'wtbx_derived', 'tif')
        stream_output_ext = "parquet" if configs.streams_as_parquet else "gpkg"
        stream_info_ext = "parquet" if configs.use_parquet else "csv"
        self.new_StrmShp_matched = self.file_path(self.strm_folder, 'matched', stream_output_ext)
        self.stream_info_file = self.file_path(self.FLDPLN_Folder, 'stream_info', stream_info_ext)
        self.new_stream_raster = self.file_path(self.strm_folder, 'matched', 'tif')
        self.lake_raster = self.file_path(self.strm_folder, 'lakes', 'tif')

        # currently the land file will be the same regardless of the streamflow source
        self.LAND_File = self.file_path(self.land_folder, 'LAND_Raster', 'tif')

        #Datasets to be Created
        streamflow_source: StreamflowSource = configs.streamflow_source
        self.DEM_StrmShp = self.file_path(self.strm_folder, 'StrmShp', stream_output_ext, 'source',
                                          legacy=f"{streamflow_source}_{self.FileName}_StrmShp.{stream_output_ext}")
        self.DEM_Reanalsyis_FlowFile = self.file_path(self.FLOW_Folder, 'Reanalysis', 'csv', 'source',
                                                      legacy=f"{streamflow_source}_{self.FileName}_Reanalysis.csv")
        self.COMID_Q_File = self.file_path(self.FLOW_Folder, '2yr_flow_initial', 'csv')

        # isolating the NWM or GEOGLOWS text in the streamflow_source variable
        strm_source = self.strm_source
        # these will only vary based upon if they are NWM or GEOGLOWS
        self.ARC_FileName_Bathy = self.file_path(self.ARC_Folder, 'ARC_Input_Bathy', self.config_end, 'source',
                                                 legacy=f"{strm_source}_ARC_Input_{self.FileName}_Bathy.{self.config_end}")
        self.ARC_FileName_for_DEM_Cleaner = self.file_path(self.ARC_Folder, 'ARC_Input_InitialFlood', 'txt', 'source',
                                                           legacy=f"{strm_source}_ARC_Input_{self.FileName}_InitialFlood.txt")
        if configs.burn_streams:
            self.DEM_File_Clean = self.file_path(dem_updated, 'fixed_Clean', 'tif')
        else:
            self.DEM_File_Clean = self.file_path(dem_updated, 'Clean', 'tif')
        self.STRM_File = self.file_path(self.strm_folder, 'STRM_Raster', 'tif', 'source')
        self.STRM_File_Clean = self.file_path(self.strm_folder, 'STRM_Raster_Clean', 'tif', 'source')

        vdt_ext = configs.vdt_file_extension
        self.VDT_File_Initial = self.file_path(self.VDT_Folder, 'VDT_Database_Initial', vdt_ext, 'source')
        self.VDT_File_Bathy = self.file_path(self.VDT_Folder, 'VDT_Database_Bathy', vdt_ext, 'source')

        self.AP_File = self.file_path(self.VDT_Folder, 'AP_Database_Bathy', vdt_ext, 'source')

        self.Curve_File = self.file_path(self.VDT_Folder, 'CurveFile', 'csv', 'source')
        self.Curve_File_Initial = self.file_path(self.VDT_Folder, 'CurveFile_Initial', 'csv', 'source')
        self.Curve_File_Bathy = self.file_path(self.VDT_Folder, 'CurveFile_Bathy', 'csv', 'source')

        self.Cross_Section_File = self.file_path(self.VDT_Folder, 'XS', 'txt', 'source')
        self.Representative_Cross_Section_File = self.file_path(self.VDT_Folder, 'Representative_XS',
                                                                'parquet' if vdt_ext == 'parquet' else 'csv', 'source')

        # self.LU_and_Streams_Water_Map = self.flood_folder / f"{strm_source}_{self.FileName}_ARC_Flood_Initial.tif"
        self.bathy_water_mask = self.file_path(self.bathy_file_folder, 'water_mask', 'tif', 'source')
        self.DepthMapFile = self.file_path(self.flood_folder, 'ARC_Depth', 'tif', 'source')
        self.ARC_BathyFile = self.file_path(self.bathy_file_folder, 'ARC_Bathy', 'tif', 'source')
        self.FS_BathyFile = self.file_path(self.bathy_file_folder, 'FS_Bathy', 'tif', 'source')

        self.floodmap_id = configs.floodmap_identifier
        if self.floodmap_id:
            self.floodmap_id = f"_{self.floodmap_id}"
        else:
            self.floodmap_id = ''

        flood_map = f"{strm_source}_{self.FileName}_ARC_Flood{self.floodmap_id}"
        self.FloodMapFile = self.file_path(self.flood_folder, 'ARC_Flood', 'tif', 'source', suffix=self.floodmap_id)
        self.FloodMapFile_Initial = self.file_path(self.flood_folder, 'ARC_Flood_Initial', 'tif', 'source', suffix=self.floodmap_id,
                                                   legacy=f"{flood_map}_Initial.tif")
        self.FloodMapFile_Initial_SHP = self.file_path(self.flood_folder, 'ARC_Flood_Initial', 'shp', 'source', suffix=self.floodmap_id,
                                                       legacy=f"{flood_map}_Initial.shp")
        self.FloodMapFile_Bathy = self.file_path(self.flood_folder, 'ARC_Flood_Bathy', 'tif', 'source', suffix=self.floodmap_id,
                                                 legacy=f"{flood_map}_Bathy.tif")
        self.FloodMapFile_Bathy_SHP = self.file_path(self.flood_folder, 'ARC_Flood_Bathy', 'shp', 'source', suffix=self.floodmap_id,
                                                     legacy=f"{flood_map}_Bathy.shp")

        # these variables will have the full specifics of the streamflow source 
        self.ARC_FileName_FloodForecast = self.file_path(self.ARC_Folder, 'ARC_Input_FloodForecast', 'txt', 'source',
                                                         legacy=f"{strm_source}_ARC_Input_{self.FileName}_FloodForecast.txt")
        self.FloodDepthFile = self.file_path(self.flood_folder, 'ARC_FloodDepth', 'tif', 'source', suffix=self.floodmap_id)
        self.FloodWSEFile = self.file_path(self.flood_folder, 'ARC_FloodWSE', 'tif', 'source', suffix=self.floodmap_id)
        self.FloodVELFile = self.file_path(self.flood_folder, 'ARC_FloodVEL', 'tif', 'source', suffix=self.floodmap_id)

        if configs.mapper == Mapper.CURVE2FLOOD_FLDPLNPY:
            self.setup_fldpln_files()

        if configs.file_names or configs.short_file_names:
            self.check_file_paths_are_unique()
        if configs.move_stream_network_to_thalweg:
            self.check_whitebox_paths()

    def check_file_paths_are_unique(self):
        """Raise if file_names or short_file_names gives two of this workspace's files the same path."""
        owners = {}
        for attr, path in vars(self).items():
            if not isinstance(path, Path) or attr == 'output_dir' or attr.lower().endswith('_folder'):
                continue
            other = owners.setdefault(path, attr)
            if other != attr:
                raise ValueError(f"Watershed '{self.watershed}': 'file_names' gives {other} and {attr} the same path, {path}.")

    def check_whitebox_paths(self):
        """
        Raise if WhiteboxTools, which moving the stream network to the thalweg runs, would be given
        a path it cannot open. Checked here so a bad output_dir, name or DEM fails before any processing.
        """
        paths = [self.fixed_dem, self.filled_dem, self.flowdir, self.flowacc, self.whitebox_stream_raster, self.new_StrmShp]
        if not self.configs.burn_streams:
            # Without burning, whitebox reads the assigned DEM instead of the fixed DEM
            paths.append(self.assigned_dem)

        for path in paths:
            path = Path(path).absolute()
            found = [char for char in WHITEBOX_UNSAFE_PATH_CHARACTERS if char in str(path)]
            if found:
                raise ValueError(
                    f"Watershed '{self.watershed}': WhiteboxTools cannot open {path} because it contains "
                    f"{' and '.join(repr(char) for char in found)}. Remove "
                    f"{', '.join(repr(char) for char in WHITEBOX_UNSAFE_PATH_CHARACTERS)} from 'output_dir', "
                    f"'name' and the DEM path, or set 'move_stream_network_to_thalweg' to False."
                )

    @property
    def model_StrmShp(self) -> Path:
        """
        The stream network the models run on: the network moved to the thalweg when moving streams,
        otherwise the source network. Steps after the move should use this rather than
        DEM_StrmShp, which minimize_output_files deletes once the moved network is written.
        """
        if self.configs.move_stream_network_to_thalweg:
            return self.new_StrmShp_matched
        return self.DEM_StrmShp

    @property
    def config_end(self) -> str:
        """Return the ARC/mapper config extension for this workspace."""
        return 'yaml' if self.configs.use_yaml else 'txt'

    def setup_dem_name(self, dem: Path | None):
        """Resolve the DEM's name and the file name stem used for its outputs."""
        if dem:
            self.dem_name = Path(dem).stem
        elif self.configs.bbox:
            # the 5th decimal place corresponds to ~1m resolution, which is more precise than we need for naming purposes, so we can round to 5 decimal places for cleaner file names
            self.dem_name = f"dem_{self.configs.bbox[0]:.5f}_{self.configs.bbox[1]:.5f}_{self.configs.bbox[2]:.5f}_{self.configs.bbox[3]:.5f}"
        else:
            raise ValueError("Must specify either a DEM or a bounding box in the configs.")
        self.FileName = self.dem_name

        if self.configs.buffer:
            if not self.configs.source_dems:
                raise ValueError(
                    "Buffering requested but no source DEMs assigned."
                )
            self.FileName += '_buffered'

    def folder_path(self, folder: str) -> Path:
        """
        Where ``folder``, one of WORKSPACE_FOLDERS, goes: its template in configs.folder_paths, else
        the 'default' template there, else '{name}/{folder}'. Templates can use {output_dir}, {name}
        (the watershed name), {dem} (the DEM's name) and {folder}; relative paths are relative to output_dir.
        """
        folder_paths = self.configs.folder_paths or {}
        template = str(folder_paths.get(folder, folder_paths.get('default', DEFAULT_FOLDER_PATH)))
        try:
            path = Path(template.format(output_dir=Path(self.configs.output_dir).absolute(), name=self.watershed,
                                        dem=self.dem_name, folder=folder))
        except (KeyError, IndexError, ValueError) as exc:
            raise ValueError(
                f"Watershed '{self.watershed}': 'folder_paths' template {template!r} for {folder} is invalid ({exc!r}). "
                "Templates can use {output_dir}, {name}, {dem} and {folder}."
            ) from exc
        return path if path.is_absolute() else Path(self.configs.output_dir) / path

    def file_path(self, folder: Path, blurb: str, ext: str, prefix: str | None = 'dem', suffix: str = '',
                  legacy: str | None = None) -> Path:
        """
        A file in ``folder`` named '<prefix><blurb><suffix>.<ext>'. ``blurb``, one of FILE_BLURBS, says what
        the file is and is what configs.file_names renames; ``suffix`` is the part that varies between files
        of that kind, like '_rp100'. ``prefix`` is 'source' for '<source>_<FileName>_', 'dem' for
        '<FileName>_' or None, and configs.short_file_names drops it. ``legacy`` is the name the file had
        before these options, for files that don't follow that pattern, and is kept unless they change it.
        An empty name in configs.file_names drops the blurb, naming the file by its suffix alone, like
        'rp100', or by its blurb if it has no suffix.
        """
        renamed = (self.configs.file_names or {}).get(blurb)
        if legacy is not None and renamed is None and not self.configs.short_file_names:
            return folder / legacy
        name = suffix.lstrip('_') if renamed == '' and suffix.strip('_') else f"{renamed or blurb}{suffix}"
        if not self.configs.short_file_names:
            name = {'source': f"{self.strm_source}_{self.FileName}_", 'dem': f"{self.FileName}_", None: ''}[prefix] + name
        return folder / f"{name}.{ext}"

    def return_period_flow_file(self, return_period: int) -> Path:
        """The flow file for one return period."""
        return self.file_path(self.FLOW_Folder, 'rp', 'csv', suffix=str(return_period))

    def setup_dem(self, dem: Path | None):
        """Resolve the assigned working DEM path."""
        self.original_dem = dem
        if dem and not self.configs.bbox and not self.configs.buffer:
            self.assigned_dem = Path(dem)
        else:
            ext = 'vrt' if self.configs.use_vrt else 'tif'
            self.assigned_dem = self.file_path(self.DEM_folder, 'dem', ext, legacy=f"{self.FileName}.{ext}")

    def setup_fldpln_files(self):
        """Add FLDPLN-specific working files to the workspace."""
        self.fldpln_library = self.file_path(self.FLDPLN_Folder, 'fldpln_library', 'parquet', None)

    def __repr__(self):
        return f"Workspace(watershed={self.watershed}, output_dir={self.output_dir})"
