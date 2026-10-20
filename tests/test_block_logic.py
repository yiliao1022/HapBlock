#!/usr/bin/env python3
"""Regression tests for block merging and gap handling."""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from hapblock_v17lib.blocks import WindowCall, label_similarity, merge_window_calls
from hapblock_v17lib.cluster import ClusterResult


def cluster(labels):
    numeric = np.asarray(labels, dtype=np.int32)
    text = [f"Hap{x:02d}" for x in labels]
    return ClusterResult(text, numeric, numeric.copy(), {}, len(set(labels)), len(labels), 1.0,
                         0.0, 0.0, 10.0, 1e13, 1.0, "PASS")


def call(start, labels=None, status="PASS"):
    return WindowCall("Chr01", start, start + 100_000, start, start + 100_000,
                      100_000, 100_000, 0.05, 20, 20, status,
                      cluster(labels) if labels is not None else None)


stable = [1, 1, 2, 2]
changed = [1, 2, 1, 2]

assert label_similarity(stable, [2, 2, 1, 1]) == 1.0
assert label_similarity(stable, changed) < 0.0

# One low-information 100-kb window should remain an annotated internal gap.
blocks = merge_window_calls(
    [call(0, stable), call(100_000, stable), call(200_000, None, "LOW_INFORMATION"),
     call(300_000, stable), call(400_000, stable)],
    ari_threshold=0.85, max_gap_bp=150_000, max_gap_fraction=0.50, min_block_bp=100_000,
)
assert len(blocks) == 1, blocks
assert blocks[0].start0 == 0 and blocks[0].end0 == 500_000
assert blocks[0].gap_intervals == [(200_000, 300_000, "BRIDGED_LOW_INFORMATION_OR_DISCORDANT_WINDOWS")]

# A sustained, incompatible membership change must create a boundary.
blocks = merge_window_calls(
    [call(0, stable), call(100_000, stable), call(200_000, changed), call(300_000, changed)],
    ari_threshold=0.85, max_gap_bp=0, max_gap_fraction=0.0, min_block_bp=100_000,
)
assert len(blocks) == 2, blocks
assert [(x.start0, x.end0) for x in blocks] == [(0, 200_000), (200_000, 400_000)]

print("PASS block merging and short-gap logic")
