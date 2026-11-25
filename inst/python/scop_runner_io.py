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



def read_staged_tsv(path: str | Path):
    import pandas as pd

    # R stages missing fields as empty strings; literal IDs such as "NULL"
    # and "001" must survive pandas' default NA and type inference.
    return pd.read_csv(
        path, sep="\t", dtype=str, keep_default_na=False, skip_blank_lines=False
    )


def validate_ids(values: list[str], label: str) -> list[str]:
    if not values or any(not value for value in values) or len(values) != len(set(values)):
        raise ValueError(f"{label} must contain unique non-empty IDs")
    return values


def read_ids(path: Path, label: str) -> list[str]:
    frame = read_staged_tsv(path)
    if frame.shape[1] != 1:
        raise ValueError(f"{label} must contain exactly one column")
    return validate_ids(frame.iloc[:, 0].tolist(), label)



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
