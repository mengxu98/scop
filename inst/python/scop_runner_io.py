from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


def read_json(path: Path, object_label: str | None = None) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if object_label is not None and not isinstance(value, dict):
        raise ValueError(f"{object_label} configuration must be a JSON object")
    return value


def write_manifest(value: dict, path: Path) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, allow_nan=False)
    temporary.replace(path)



def read_ids(path: Path, label: str) -> list[str]:
    import pandas as pd

    frame = pd.read_csv(path, sep="\t", dtype=str)
    if frame.shape[1] != 1:
        raise ValueError(f"{label} must contain exactly one column")
    values = frame.iloc[:, 0].astype(str).tolist()
    if not values or len(values) != len(set(values)) or any(not value for value in values):
        raise ValueError(f"{label} must contain unique non-empty IDs")
    return values



def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()



def configuration_path(description: str | None = None) -> Path:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--config", required=True, type=Path)
    return parser.parse_args().config
