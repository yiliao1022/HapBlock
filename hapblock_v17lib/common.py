from __future__ import annotations

import csv
import gzip
import json
import math
import re
from pathlib import Path
from typing import Iterable, Sequence


def open_text(path, mode="rt"):
    path = str(path)
    return gzip.open(path, mode, encoding="utf-8") if path.endswith(".gz") else open(path, mode, encoding="utf-8")


def read_list(path) -> list[str]:
    result = []
    seen = set()
    if not path:
        return result
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        value = line.strip().split()[0] if line.strip() and not line.lstrip().startswith("#") else ""
        if value and value not in seen:
            result.append(value)
            seen.add(value)
    return result


def write_tsv(path, fields: Sequence[str], rows: Iterable[dict]):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def chrom_key(chrom: str):
    match = re.search(r"(\d+)$", chrom)
    return (chrom[:match.start()], int(match.group(1))) if match else (chrom, 10**12)


def fi(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def ff(value, default=math.nan):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_int_csv(text: str) -> list[int]:
    values = [int(x.strip()) for x in text.split(",") if x.strip()]
    if not values or any(x <= 0 for x in values):
        raise ValueError("Window/step values must be positive integers")
    return values


def parse_float_csv(text: str) -> list[float]:
    values = [float(x.strip()) for x in text.split(",") if x.strip()]
    if not values:
        raise ValueError("Expected at least one numeric value")
    return values


def stable_subsample_indices(n: int, cap: int):
    import numpy as np
    if cap <= 0 or n <= cap:
        return np.arange(n, dtype=np.int64)
    return np.unique(np.linspace(0, n - 1, cap, dtype=np.int64))
