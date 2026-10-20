#!/usr/bin/env python3
"""Cache reuse must be deterministic and reject mismatched sample lists."""
import tempfile
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from hapblock_v17lib.cache import build_cache


testdir = ROOT / "tests"
with tempfile.TemporaryDirectory() as directory:
    directory = Path(directory)
    cache = directory / "cache"
    first = build_cache(testdir / "synthetic.vcf", testdir / "cultivar.list", testdir / "wild.list", cache)
    second = build_cache(testdir / "synthetic.vcf", testdir / "cultivar.list", testdir / "wild.list", cache, resume=True)
    assert first == second
    assert len(first["sample_qc"]) == 44
    assert all(0.0 <= row["call_rate"] <= 1.0 for row in first["sample_qc"])

    changed = directory / "changed.list"
    changed.write_text("C00\n", encoding="utf-8")
    try:
        build_cache(testdir / "synthetic.vcf", changed, testdir / "wild.list", cache, resume=True)
    except SystemExit as error:
        assert "sample lists differ" in str(error)
    else:
        raise AssertionError("mismatched sample list unexpectedly reused cache")

print("PASS cache fingerprint and sample QC")
