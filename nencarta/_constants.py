from pathlib import Path

THIS_DIR = Path(__file__).parent
ESA_TILES_FILE = THIS_DIR / 'data' / 'esa_tiles.gpkg'
GEOGLOWS_RETURN_PERIODS_URL = 's3://geoglows-v2/retrospective/return-periods.zarr'
GEOGLOWS_FDC_URL = 's3://geoglows-v2/retrospective/fdc.zarr'
GEOGLOWS_DAILY_URL = 's3://geoglows-v2/retrospective/daily.zarr'
GEOGLOWS_FORECAST_PREFIX_URL = 's3://geoglows-v2-forecasts/'
NWM_RP_URL = 'https://nwm-api.ciroh.org/return-period'
CACHED_ESA_GRID = THIS_DIR / '.esa_worldcover_grid.parquet'
FLOODSPREADER_PATH = THIS_DIR / "floodspreader.py"