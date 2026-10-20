#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from hapblock_v17lib.cluster import DistanceResult, association_status, complete_linkage_clusters


def main():
    # A-B and B-C pass, A-C fails. Complete linkage must not put all three together.
    distance = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]], dtype=float)
    counts = np.ones((3, 3), dtype=np.int64) * 10
    result = complete_linkage_clusters(DistanceResult(distance, counts, counts, 10, 10, 1000000),
                                       ["A", "B", "C"], 1.0, 0.0, 1)
    assert len(set(result.raw_labels.tolist())) == 2, result.raw_labels
    labels = ["Hap01"] * 120 + ["Hap02"] * 80 + ["Hap03"] * 52 + ["Hap04"] * 44 + ["Minor01"] * 104
    status = association_status(labels, 400, 0.10, 0.50, 2)
    assert status["association_ready"]
    assert status["n_major_haplotypes"] == 4
    assert abs(status["major_assigned_fraction"] - 0.74) < 1e-12
    print("PASS cluster and association rules")


if __name__ == "__main__":
    main()
