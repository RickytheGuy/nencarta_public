import threading
import os
import pickle
from typing import Any
from pathlib import Path

_HAS_LMDB = False
try:
    # If lmdb is avaialable, we can significantly speed up raster metadata loading by caching it in an lmdb database.
    # This is especially good if there are thousands of rasters to process
    import lmdb
    _HAS_LMDB = True
except (ImportError, ModuleNotFoundError):
    pass


class LMDBCache:
    def __init__(self, cache_name: str):
        self.can_cache = _HAS_LMDB
        if self.can_cache:
            self._setup_cache(cache_name)
        

    def _setup_cache(self, cache_name: str) -> None:
        self._CACHE_LOCK = threading.RLock()

        self._CACHE_DIR = Path.home() / ".cache"
        self._CACHE_DIR.mkdir(parents=True, exist_ok=True)

        LMDB_PATH = self._CACHE_DIR / f"{cache_name}.lmdb"
        self._LMDB_ENV = lmdb.open(
            str(LMDB_PATH),
            map_size=1024**3,  # 1 GB
            subdir=False,
            lock=True,
            sync=False,
            metasync=False,
            readahead=False,
            writemap=True,
            max_readers=512,
        )

    def _save_cached_metadata(
        self,
        filepath: str,
        timestamp: float,
        metadata: dict[str, Any],
    ) -> None:

        key = self._cache_key(filepath)

        value = pickle.dumps(
            {
                "timestamp": timestamp,
                "metadata": metadata,
            },
            protocol=pickle.HIGHEST_PROTOCOL,
        )

        with self._LMDB_ENV.begin(write=True) as txn:
            txn.put(key, value)

    def _load_cached_metadata(
        self,
        filepath: str,
        timestamp: float,
    ) -> dict[str, Any] | None:

        key = self._cache_key(filepath)

        with self._LMDB_ENV.begin(write=False) as txn:
            value = txn.get(key)

        if value is None:
            return None

        try:
            entry = pickle.loads(value)
        except Exception:
            return None

        if entry["timestamp"] != timestamp:
            return None

        return entry["metadata"]
    
    def _get_file_timestamp(self, filepath: str) -> float:
        return os.path.getmtime(filepath)
    
    def _cache_key(self, filepath: str) -> bytes:
        return os.path.abspath(filepath).encode()
