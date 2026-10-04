import hashlib
import json
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
IGNORE_DIR = REPO_ROOT / "ignore"
BASELINE_MANIFEST = REPO_ROOT / "tests" / "fixtures" / "ignore_output_manifest.json"

SOUTH_PLATTE_INPUTS = [
    Path(r"C:\Users\lrr43\Documents\masters\fim_sites\South_Platte_River_at_Fort_Morgan_(2017)\streams.parquet"),
    Path(r"C:\Users\lrr43\Documents\masters\fim_sites\South_Platte_River_at_Fort_Morgan_(2017)\DEMs"),
    Path(r"C:\Users\lrr43\Documents\masters\fim_sites\South_Platte_River_at_Fort_Morgan_(2017)\flow_files\stage=20.csv"),
]

IGNORED_OUTPUT_FOLDERS = [
    "forecast_test_original",
    "forecast_test_new",
    "clean_dem_original",
    "clean_dem_new",
    "no_bathy_original",
    "no_bathy_new",
]


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_manifest(root: Path) -> dict:
    files = {}
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        relative_path = path.relative_to(root).as_posix()
        files[relative_path] = {
            "size": path.stat().st_size,
            "sha256": file_sha256(path),
        }
    return files


@pytest.mark.parametrize("path", SOUTH_PLATTE_INPUTS)
def test_south_platte_inputs_used_by_ignore_calls_still_exist(path):
    assert path.exists()


@pytest.mark.parametrize("folder_name", IGNORED_OUTPUT_FOLDERS)
def test_cached_output_folders_exist(folder_name):
    assert (IGNORE_DIR / folder_name).is_dir()


def test_cached_ignore_outputs_match_manifest():
    if not BASELINE_MANIFEST.exists():
        pytest.skip("Run tests/write_ignore_output_manifest.py to create the local baseline.")

    expected = json.loads(BASELINE_MANIFEST.read_text())
    actual = {
        folder_name: build_manifest(IGNORE_DIR / folder_name)
        for folder_name in IGNORED_OUTPUT_FOLDERS
    }

    assert actual == expected
