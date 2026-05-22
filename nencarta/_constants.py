import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ESA_TILES_FILE = os.path.join(THIS_DIR, 'data', 'esa_tiles.gpkg')
GEOGLOWS_RETURN_PERIODS_URL = 's3://geoglows-v2/retrospective/return-periods.zarr'
GEOGLOWS_FDC_URL = 's3://geoglows-v2/retrospective/fdc.zarr'
GEOGLOWS_DAILY_URL = 's3://geoglows-v2/retrospective/daily.zarr'
NWM_RP_URL = 'https://nwm-api.ciroh.org/return-period'
CACHED_ESA_GRID = os.path.join(os.path.dirname(__file__), '.esa_worldcover_grid.parquet')
FLOODSPREADER_PATH = os.path.join(THIS_DIR, "floodspreader.py")